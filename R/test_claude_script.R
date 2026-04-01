# =============================================================================
# Joint Species Distribution Model (JSDM) for Tadpole Gut Microbiome Data
# Package: Hmsc
# Data: ASV read counts (samples x ASVs), metadata, per-area environmental data
# =============================================================================

# =============================================================================
# SECTION 1: Install and load required packages
# =============================================================================

# Run this block once to install packages if you don't have them yet.
# After the first install, you can comment these lines out.

library(Hmsc)
library(tidyverse)
library(corrplot)


# =============================================================================
# SECTION 2: Load your data files
# =============================================================================

# --- 2a. Load the combined ASV count table ---
# Replace the filename with your actual file path/name.
# This file should have samples as rows and ASVs as columns,
# with the first column being sample IDs.

asv_raw <- read.csv("data/raw/otu_table.csv",      
                    row.names = 1,            # sets the first column as row names (sample IDs)
                    check.names = FALSE)      # preserves ASV names exactly as-is

# Quick check: dimensions should be (number of samples) x (number of ASVs)
cat("ASV table dimensions (samples x ASVs):", dim(asv_raw), "\n")


# --- 2b. Load the metadata file ---
# This file should have one row per sample, with columns for sample ID,
# area/pond info, and any other sample-level variables.
# Replace the filename with your actual file path/name.

metadata <- read.csv("data/raw/DCNF_metadata.csv",         
                     stringsAsFactors = FALSE)

# Tell R which column in metadata holds the sample IDs.
# Replace "SampleID" with your actual column name.
sample_id_col <- "Sample_name"

# Quick check
cat("Metadata dimensions:", dim(metadata), "\n")
cat("Metadata columns:", colnames(metadata), "\n")


# --- 2c. Load the per-area environmental data ---
# You have one file per area. Each file has samples as rows and
# environmental variables as columns (temp, pH, tadpole species counts, etc.).
# List your three files here with their actual filenames.

#env_area1 <- read.csv("env_area1.csv",       
                      #stringsAsFactors = FALSE)
#env_area2 <- read.csv("env_area2.csv",       
                      #stringsAsFactors = FALSE)
#env_area3 <- read.csv("env_area3.csv",       
                      #stringsAsFactors = FALSE)

# Add an "Area" label column to each, matching the area name in your metadata.
# Replace "Area1", "Area2", "Area3" with whatever your areas are actually called.
#env_area1$Area <- "Area1"                    # <-- CHANGE to your area name
#env_area2$Area <- "Area2"                    # <-- CHANGE to your area name
#env_area3$Area <- "Area3"                    # <-- CHANGE to your area name

# Combine all three environmental files into one
#env_all <- bind_rows(env_area1, env_area2, env_area3)

#cat("Combined environmental data dimensions:", dim(env_all), "\n")
#cat("Environmental columns:", colnames(env_all), "\n")


# =============================================================================
# SECTION 3: Fix the hyphen mismatch between ASV table and metadata
# =============================================================================
# In your ASV file, some sample IDs had hyphens removed (e.g. "SA001" instead
# of "SA-001"). We standardise both files by removing all hyphens from IDs
# before joining, so everything lines up correctly.

# Remove hyphens from the ASV table row names
rownames(asv_raw) <- gsub("-", "", rownames(asv_raw))

# Remove hyphens from the metadata sample ID column
metadata[[sample_id_col]] <- gsub("-", "", metadata[[sample_id_col]])

# Also remove hyphens from the sample ID column in the environmental data.
# Replace "SampleID" below with the actual sample ID column name in your env files.
#env_id_col <- "SampleID"                     # <-- CHANGE if different in env files
#env_all[[env_id_col]] <- gsub("-", "", env_all[[env_id_col]])

cat("Sample IDs in ASV table (first 5):", head(rownames(asv_raw), 5), "\n")
cat("Sample IDs in metadata (first 5):", head(metadata[[sample_id_col]], 5), "\n")


# =============================================================================
# SECTION 4: Choose your area of interest and filter data
# =============================================================================
# Set the area you want to model. This must match the "Area" label you used
# in Section 2c above.

#target_area <- "Area1"                       # <-- CHANGE to your area of interest

# Filter metadata to only samples from that area.
# Replace "Area" below with the actual area column name in your metadata.
#area_col_metadata <- "Area"                  # <-- CHANGE to your area column name in metadata

#meta_area <- metadata %>%
  #filter(.data[[area_col_metadata]] == target_area)

#cat("Number of samples in", target_area, ":", nrow(meta_area), "\n")

# Filter the ASV table to only the samples in that area
# (keeping only rows whose IDs appear in the filtered metadata)
samples_in_area <- metadata[[sample_id_col]]
asv_dcnf <- asv_raw[rownames(asv_raw) %in% samples_in_area, ]

cat("ASV table filtered to area:", dim(asv_area), "\n")

# Filter the environmental data to only samples in that area
#env_area <- env_all %>%
  #filter(Area == target_area)

#cat("Environmental data filtered to area:", dim(env_area), "\n")


# =============================================================================
# SECTION 5: Select the top 50 ASVs by total read count in this area
# =============================================================================

# Sum reads per ASV across all samples in this area
asv_totals <- colSums(asv_dcnf)

# Sort descending and take the top 50
top50_names <- names(sort(asv_totals, decreasing = TRUE))[1:50]

cat("Top 50 ASVs selected. Most abundant:", top50_names[1:5], "\n")

# Subset the ASV table to just those 50 ASVs
asv_top50 <- asv_dcnf[, top50_names]


# =============================================================================
# SECTION 6: Align all data frames so rows match in the same order
# =============================================================================
# All three data objects (ASV counts, metadata, environmental data) must have
# exactly the same samples in exactly the same row order for Hmsc to work.

# Use the ASV table row order as the reference
sample_order <- rownames(asv_top50)

# Reorder metadata rows to match
meta_aligned <- metadata %>%
  column_to_rownames(sample_id_col) %>%   # make sample IDs the row names
  slice(match(sample_order, rownames(.)))

# Reorder environmental data rows to match
#env_aligned <- env_area %>%
  #column_to_rownames(env_id_col) %>%
  #slice(match(sample_order, rownames(.)))

# Remove the "Area" column from env data (it's not a numeric predictor)
#env_aligned <- env_aligned %>% select(-Area)

# Confirm alignment
cat("Rows aligned — ASV:", nrow(asv_top50),
    "| Metadata:", nrow(meta_aligned), "\n")

# =============================================================================
# SECTION 7: Prepare HMSC model inputs
# =============================================================================

# --- 7a. Y matrix: the response variable (ASV counts) ---
# Hmsc expects a numeric matrix. Rows = samples, Columns = ASVs.
Y <- as.matrix(asv_top50)

# --- 7b. XData: the environmental predictors (fixed effects) ---
# This is your environmental data frame. All columns should be numeric.
# List the column names you want to use as predictors.
# Replace these with your actual environmental column names.
predictor_cols <- c("watertemp", "conductivity", "chlorophyll",
                    "DO", "pH", "Gosner_stage")

XData <- meta_aligned %>%
  select(all_of(predictor_cols)) %>%
  mutate(across(everything(), as.numeric))  # ensure all are numeric

# Remove rows with NAs and track which samples were kept
complete_rows <- complete.cases(XData)
XData        <- XData[complete_rows, ]
Y            <- Y[complete_rows, ]

studyDesign <- data.frame(
  pond = as.factor(meta_aligned[[pond_col]][complete_rows])
)


# Scale predictors to mean=0, sd=1 (strongly recommended for HMSC)
XData <- as.data.frame(scale(XData))

# --- 7c. XFormula: tells Hmsc which predictors to use ---
# This uses all columns in XData. You can modify this formula to add
# interactions (e.g. temperature * pH) if needed.
XFormula <- as.formula(paste("~", paste(colnames(XData), collapse = " + ")))
cat("Model formula:", deparse(XFormula), "\n")

# --- 7d. studyDesign: the random effects structure ---
# Here we use "pond" as a random effect (samples from the same pond
# are not independent). Replace "Pond" with your actual pond column name.
# If you don't have a pond column, you can use sample ID itself.
pond_col <- "Pond"                      

studyDesign <- data.frame(
  pond = as.factor(meta_aligned[[pond_col]])
)

# Create the random level object for pond
rL_pond <- HmscRandomLevel(units = levels(studyDesign$pond))


# =============================================================================
# SECTION 8: Define and fit the HMSC model
# =============================================================================

# Define the model
# distr = "lognormal poisson" is appropriate for raw integer read counts
# (it accounts for the overdispersion typical in microbiome data).
# Alternatives: "poisson" (simpler), "normal" (for log-transformed data).

model <- Hmsc(
  Y           = Y,
  XData       = XData,
  XFormula    = XFormula,
  studyDesign = studyDesign,
  ranLevels   = list(pond = rL_pond),
  distr       = "lognormal poisson"
)

cat("Hmsc model object created successfully.\n")

# --- Fit the model using MCMC sampling ---
# These settings are for a QUICK TEST RUN to make sure everything works.
# For a real analysis you should use much higher values (see commented block below).

nChains   <- 2       # number of MCMC chains (use 4 for final run)
thin      <- 1       # thinning interval (use 10-100 for final run)
samples   <- 100     # posterior samples per chain (use 1000+ for final run)
transient <- 50      # burn-in samples to discard (use 500+ for final run)

set.seed(42)         # for reproducibility

m_fitted <- sampleMcmc(
  model,
  nChains   = nChains,
  thin      = thin,
  samples   = samples,
  transient = transient,
  verbose   = 50      # prints progress every 50 iterations
)

cat("Model fitting complete!\n")

# --- Settings for a FINAL / PUBLICATION-QUALITY run ---
# Uncomment and use these instead of the test settings above when ready.
# Expect this to take several hours depending on your machine.
#
# nChains   <- 4
# thin      <- 10
# samples   <- 1000
# transient <- 500


# =============================================================================
# SECTION 9: Check MCMC convergence
# =============================================================================
# Rhat values close to 1.0 (ideally < 1.1) indicate good convergence.
# Effective sample sizes (ESS) should be reasonably large (>100).

mpost <- convertToCodaObject(m_fitted)

# Check convergence of the beta parameters (species-environment relationships)
beta_psrf <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
cat("Beta Rhat values (should be close to 1):\n")
print(summary(beta_psrf[, 1]))   # shows min, mean, max Rhat


# =============================================================================
# SECTION 10: Examine and visualise results
# =============================================================================

# --- 10a. Variance partitioning ---
# Shows how much of the variation in each ASV is explained by each
# environmental predictor vs. the random effect (pond).

VP <- computeVariancePartitioning(m_fitted)

# Plot variance partitioning
plotVariancePartitioning(m_fitted, VP = VP,
                         main = paste("Variance Partitioning —", "DCNF"))


# --- 10b. Species-environment relationships (beta parameters) ---
# Shows which environmental variables are associated with each ASV.
# Red = positive association, blue = negative association.
# Only relationships with >95% posterior support are shown by default.

postBeta <- getPostEstimate(m_fitted, parName = "Beta")

plotBeta(m_fitted,
         post  = postBeta,
         param = "Support",          # shows posterior support (confidence)
         supportLevel = 0.95,        # only show effects with 95%+ support
         main  = paste("Species-Environment Relationships —", DCNF))


# --- 10c. Species associations (omega / residual co-occurrence) ---
# Shows which pairs of ASVs tend to co-occur (positive, red) or
# avoid each other (negative, blue) after accounting for environment.

OmegaCor <- computeAssociations(m_fitted)

# Plot the association matrix for the first random level (pond)
supportLevel <- 0.95

toPlot <- ((OmegaCor[[1]]$support > supportLevel)
           + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

corrplot(toPlot,
         method  = "color",
         col     = colorRampPalette(c("blue", "white", "red"))(200),
         title   = paste("Residual ASV Associations —", target_area),
         mar     = c(0, 0, 2, 0),
         tl.cex  = 0.5,        # shrinks ASV labels so they fit
         tl.col  = "black")


# --- 10d. Model fit (R-squared per ASV) ---
# Shows how well the model predicts each ASV's abundance.

preds      <- computePredictedValues(m_fitted)
model_fit  <- evaluateModelFit(hM = m_fitted, predY = preds)

cat("Mean Tjur R2 across top 50 ASVs:", mean(model_fit$TjurR2, na.rm = TRUE), "\n")
cat("Per-ASV Tjur R2:\n")
print(data.frame(ASV = colnames(Y), TjurR2 = round(model_fit$TjurR2, 3)))


# =============================================================================
# SECTION 11: Save the fitted model
# =============================================================================
# Saving means you don't have to rerun the slow MCMC step every time.
# Load it back with: m_fitted <- readRDS("hmsc_model_area1.rds")

saveRDS(m_fitted, file = paste0("hmsc_model_", "DCNF", ".rds"))
cat("Model saved to hmsc_model_", "DCNF", ".rds\n", sep = "")
