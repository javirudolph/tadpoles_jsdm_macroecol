# Load packages
library(jSDM)
library(dplyr)
library(vegan)
library(stringr)

# Load Data
meta <- read.csv("data/raw/DCNF_metadata.csv")
otu  <- read.csv("data/raw/otu_table.csv", row.names = 1)
tax  <- read.csv("data/raw/tax_table.csv")

# Step 2: Standardize sample names (remove hyphens from metadata)
meta <- meta |> 
  mutate(sample_id_clean = gsub("-", "", Sample_name))

# Step 3: Define DCNF ponds
dcnf_ponds <- c("08", "181", "191", "291", "292", "294", 
                "361", "871", "SFA2", "SFA4")

# Step 4: Extract pond IDs by searching for pond codes in sample names
sample_names <- rownames(otu)

# Function to find which pond code is in each sample name
extract_pond <- function(sample_name, pond_list) {
  for (pond in pond_list) {
    if (grepl(pond, sample_name)) {
      return(pond)
    }
  }
  return(NA)
}

# Apply to all sample names
pond_ids <- sapply(sample_names, extract_pond, pond_list = dcnf_ponds)

# Create lookup table
sample_lookup <- data.frame(
  sample_id_otu = sample_names,
  pond_id = pond_ids,
  stringsAsFactors = FALSE
)

# Step 5: Filter to DCNF samples (non-NA pond IDs)
dcnf_sample_ids <- sample_lookup %>%
  filter(!is.na(pond_id)) %>%
  pull(sample_id_otu)

otu_dcnf <- otu[dcnf_sample_ids, ]

cat("Number of DCNF samples in OTU table:", nrow(otu_dcnf), "\n")
cat("Ponds represented:", unique(sample_lookup$pond_id[!is.na(sample_lookup$pond_id)]), "\n")

# Step 6: Match with metadata
metadata_dcnf <- meta |> 
  filter(sample_id_clean %in% dcnf_sample_ids) |> 
  arrange(match(sample_id_clean, rownames(otu_dcnf)))

all(rownames(otu_dcnf) == metadata_dcnf$sample_id_clean)

# Step 7: Add pond_id to metadata if not already there
metadata_dcnf <- metadata_dcnf %>%
  left_join(sample_lookup, by = c("sample_id_clean" = "sample_id_otu"))

# Verify
cat("\nSamples in OTU table:", nrow(otu_dcnf), "\n")
cat("Samples in metadata:", nrow(metadata_dcnf), "\n")
cat("Ponds in metadata:", length(unique(metadata_dcnf$pond_id)), "\n")
cat("Samples per pond:\n")
print(table(metadata_dcnf$pond_id))

# Step 8: Filter rare SVs
# Keep SVs present in at least 10% of samples
prevalence_threshold <- 0.10
min_samples <- ceiling(prevalence_threshold * nrow(otu_dcnf))

sv_prevalence <- colSums(otu_dcnf > 0)
svs_to_keep <- sv_prevalence >= min_samples

otu_filtered <- otu_dcnf[, svs_to_keep]

cat("Original SVs:", ncol(otu_dcnf), "\n")
cat("SVs after filtering:", ncol(otu_filtered), "\n")
cat("Prevalence threshold:", min_samples, "samples\n")

# Step 9: Convert to presence-absence
pa_matrix <- (otu_filtered > 0) * 1

# Step 10: Prepare environmental data for simple model
# Let's start with: Family, Type (fish), Gosner stage, and one abiotic variable (pH)
env_data <- metadata_dcnf %>%
  mutate(
    Family = as.factor(Family),
    Type = as.factor(Type),
    Gosner_stage = as.numeric(Gosner_stage),
    pH = as.numeric(pH),
    pond_id = as.factor(pond_id)
  ) %>%
  select(Family, Type, Gosner_stage, pH, pond_id)

# Check for missing values
cat("\nMissing values per variable:\n")
print(colSums(is.na(env_data)))

# Remove any samples with missing data (if needed)
complete_cases <- complete.cases(env_data)
if(sum(!complete_cases) > 0) {
  cat("\nRemoving", sum(!complete_cases), "samples with missing data\n")
  pa_matrix <- pa_matrix[complete_cases, ]
  env_data <- env_data[complete_cases, ]
}

# Step 11: Scale continuous variables
env_data_scaled <- env_data %>%
  mutate(
    Gosner_stage_scaled = as.numeric(scale(Gosner_stage)),
    pH_scaled = as.numeric(scale(pH))
  )

# Step 12: Simple model formula
# Starting with main effects only
formula_simple <- ~ Family + Type + Gosner_stage_scaled + pH_scaled

# Step 13: Fit jSDM
# Starting with 2 latent variables to capture co-occurrence patterns
cat("\nFitting jSDM model...\n")
cat("This may take several minutes...\n\n")

model_simple <- jSDM_binomial_probit(
  presence_data = pa_matrix,
  site_formula = formula_simple,
  site_data = env_data_scaled,
  n_latent = 2,  # latent variables for species co-occurrence
  burnin = 2000,
  mcmc = 2000,
  thin = 2,
  alpha_start = 0,
  beta_start = 0,
  lambda_start = 0,
  W_start = 0,
  V_alpha = 1,
  shape_Valpha = 0.5,
  rate_Valpha = 0.0005,
  mu_beta = 0,
  V_beta = 1,
  mu_lambda = 0,
  V_lambda = 1,
  seed = 123,
  verbose = 1
)

cat("\nModel fitting complete!\n")

# Step 14: Examine model results

# Summary of the model
summary(model_simple)

# Step 15: Check model structure
cat("\nModel components:\n")
print(names(model_simple))

cat("\nStructure of mcmc.sp (first element):\n")
print(str(model_simple$mcmc.sp[[1]]))

# Step 16: Extract beta coefficients properly
n_species <- length(model_simple$mcmc.sp)

# Get one element to check dimensions
first_mcmc <- as.matrix(model_simple$mcmc.sp[[1]])
predictor_names <- colnames(first_mcmc)
n_predictors <- ncol(first_mcmc)

cat("\nNumber of species:", n_species, "\n")
cat("Number of predictors:", n_predictors, "\n")
cat("Predictor names:", predictor_names, "\n")

# Create beta matrix
beta_matrix <- matrix(NA, nrow = n_species, ncol = n_predictors)

for(i in 1:n_species) {
  mcmc_mat <- as.matrix(model_simple$mcmc.sp[[i]])
  beta_matrix[i, ] <- colMeans(mcmc_mat)
}

colnames(beta_matrix) <- predictor_names
rownames(beta_matrix) <- colnames(pa_matrix)

cat("\nBeta matrix created successfully!\n")
cat("Dimensions:", dim(beta_matrix), "\n\n")

# First few rows
cat("First few rows of beta matrix:\n")
print(head(beta_matrix, 3))

# Step 17: Average effect sizes
predictor_effects <- colMeans(abs(beta_matrix))
cat("\nAverage absolute effect size by predictor:\n")
print(round(predictor_effects, 3))

# paste the most recent code output here from tadpoles II