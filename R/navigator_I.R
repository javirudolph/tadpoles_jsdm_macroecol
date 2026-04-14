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

# Step 18: Extract and visualize latent variables

cat("\nNumber of species in pa_matrix:", ncol(pa_matrix), "\n")

# Extract latent variables from MCMC chains
W1_samples <- as.matrix(model_simple$mcmc.latent$lv_1)
W2_samples <- as.matrix(model_simple$mcmc.latent$lv_2)

cat("W1 dimensions:", dim(W1_samples), "\n")
cat("W2 dimensions:", dim(W2_samples), "\n")

# Use the actual number of species that were modeled
n_species_modeled <- ncol(W1_samples)

lambda_matrix <- matrix(NA, nrow = n_species_modeled, ncol = 2)
lambda_matrix[, 1] <- colMeans(W1_samples)
lambda_matrix[, 2] <- colMeans(W2_samples)

# Need to figure out which species were actually modeled
# Check if beta_matrix has the same number of species
cat("Beta matrix rows (species):", nrow(beta_matrix), "\n")

# Use species names from beta_matrix since that's what was actually modeled
if(nrow(beta_matrix) == n_species_modeled) {
  rownames(lambda_matrix) <- rownames(beta_matrix)
} else {
  rownames(lambda_matrix) <- paste0("SV_", 1:n_species_modeled)
}

colnames(lambda_matrix) <- c("Latent_Variable_1", "Latent_Variable_2")

cat("\nLambda matrix created!\n")
cat("Dimensions:", dim(lambda_matrix), "\n")
print(head(lambda_matrix))

# Plot co-occurrence structure
par(mfrow = c(1, 1))
plot(lambda_matrix[,1], lambda_matrix[,2],
     xlab = "Latent Variable 1",
     ylab = "Latent Variable 2",
     main = "SV Co-occurrence Structure",
     pch = 16, cex = 0.7, col = rgb(0, 0, 1, 0.5))
abline(h = 0, v = 0, lty = 2, col = "gray")

# Step 19: Top SVs responding to Family
family_cols <- grep("Family", colnames(beta_matrix))
if(length(family_cols) > 0) {
  family_effects <- rowSums(abs(beta_matrix[, family_cols]))
  top_family_idx <- order(family_effects, decreasing = TRUE)[1:10]
  
  cat("\n\nTop 10 SVs responding to Family:\n")
  top_svs <- data.frame(
    SV = rownames(beta_matrix)[top_family_idx],
    Effect_size = round(family_effects[top_family_idx], 3)
  )
  print(top_svs)
}

# Step 20: Top SVs responding to Fish presence (Type)
type_cols <- grep("Type", colnames(beta_matrix))
if(length(type_cols) > 0) {
  type_effects <- rowSums(abs(beta_matrix[, type_cols]))
  top_type_idx <- order(type_effects, decreasing = TRUE)[1:10]
  
  cat("\n\nTop 10 SVs responding to Fish presence:\n")
  top_svs_type <- data.frame(
    SV = rownames(beta_matrix)[top_type_idx],
    Effect_size = round(type_effects[top_type_idx], 3)
  )
  print(top_svs_type)
}

# Expanded model formula including pond
# Note: jSDM doesn't have built-in random effects, so we'll include pond as a fixed effect
formula_expanded <- ~ Family + Type + Gosner_stage_scaled + pH_scaled + pond_id

cat("\nFitting expanded jSDM model with pond...\n")
cat("This may take several minutes...\n\n")

model_expanded <- jSDM_binomial_probit(
  presence_data = pa_matrix,
  site_formula = formula_expanded,
  site_data = env_data_scaled,
  n_latent = 2,
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

cat("\nExpanded model fitting complete!\n")

# Prepare data for visualization
library(tidyverse)

# Add pond and family info to the presence-absence matrix
pa_with_info <- pa_matrix %>%
  as.data.frame() %>%
  mutate(
    pond_id = env_data_scaled$pond_id,
    Family = env_data_scaled$Family
  )

# Calculate SV richness (number of SVs present) per sample
sv_richness <- data.frame(
  pond_id = env_data_scaled$pond_id,
  Family = env_data_scaled$Family,
  SV_richness = rowSums(pa_matrix)
)

# Bar plot: SV richness by pond
plot_pond <- ggplot(sv_richness, aes(x = pond_id, y = SV_richness, fill = pond_id)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "SV Richness by Pond",
       x = "Pond ID",
       y = "Number of SVs Present") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_pond)

# Bar plot: SV richness by family
plot_family <- ggplot(sv_richness, aes(x = Family, y = SV_richness, fill = Family)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "SV Richness by Tadpole Family",
       x = "Family",
       y = "Number of SVs Present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_family)

# Composite plot: SV richness by family within each pond
plot_combined <- ggplot(sv_richness, aes(x = Family, y = SV_richness, fill = Family)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  facet_wrap(~ pond_id, ncol = 5) +
  labs(title = "SV Richness by Family across Ponds",
       x = "Family",
       y = "Number of SVs Present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 9))

print(plot_combined)

# Summary statistics
cat("\nSV Richness Summary by Pond:\n")
print(sv_richness %>% 
        group_by(pond_id) %>% 
        summarise(
          mean_richness = round(mean(SV_richness), 1),
          sd_richness = round(sd(SV_richness), 1),
          n_samples = n()
        ))

cat("\nSV Richness Summary by Family:\n")
print(sv_richness %>% 
        group_by(Family) %>% 
        summarise(
          mean_richness = round(mean(SV_richness), 1),
          sd_richness = round(sd(SV_richness), 1),
          n_samples = n()
        ))

# Calculate prevalence of each SV by pond
sv_by_pond <- pa_with_info %>%
  pivot_longer(cols = -c(pond_id, Family), 
               names_to = "SV", 
               values_to = "presence") %>%
  group_by(pond_id, SV) %>%
  summarise(prevalence = mean(presence), .groups = "drop")

# Top 10 SVs per pond
top_svs_pond <- sv_by_pond %>%
  group_by(pond_id) %>%
  slice_max(order_by = prevalence, n = 10)

# Heatmap of top SVs by pond
ggplot(top_svs_pond, aes(x = pond_id, y = SV, fill = prevalence)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "Top 10 SVs Prevalence by Pond",
       x = "Pond ID",
       y = "SV",
       fill = "Prevalence") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6))

library(pheatmap)
library(corrplot)

# ============================================
# Panel A: Environmental Response (like your image)
# ============================================

# This shows how each SV responds to environmental variables
# We'll use the beta_matrix you already created

# Select top responding SVs for cleaner visualization
top_n_svs <- 30
top_svs_idx <- order(rowSums(abs(beta_matrix)), decreasing = TRUE)[1:top_n_svs]

# Subset to top SVs
beta_subset <- beta_matrix[top_svs_idx, ]

# Create heatmap
pheatmap(beta_subset,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-4, 4, length.out = 101),
         main = "Environmental Response",
         fontsize_row = 8,
         fontsize_col = 10)

# ============================================
# Panel B: Covariance/Correlation Matrix
# ============================================

# Calculate pairwise correlations between SVs based on their environmental responses
sv_correlations <- cor(t(beta_matrix))

# Select same top SVs for comparison
sv_cor_subset <- sv_correlations[top_svs_idx, top_svs_idx]

pheatmap(sv_cor_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "SV Response Covariance",
         fontsize_row = 8,
         fontsize_col = 8)

# ============================================
# Panel C: Co-occurrence in Data
# ============================================

# Calculate observed co-occurrence from presence-absence matrix
# Correlation between SV occurrences across samples
sv_cooccurrence <- cor(pa_matrix)

# Select top SVs
sv_cooc_subset <- sv_cooccurrence[top_svs_idx, top_svs_idx]

pheatmap(sv_cooc_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "Co-occurrence in Data",
         fontsize_row = 8,
         fontsize_col = 8)

# ============================================
# Panel D: Residual Correlation (FIXED v3)
# ============================================

# The lambda_matrix only has 76 species that were actually modeled
# We need to identify WHICH 76 from the original 405

# Check which SVs from beta_matrix correspond to lambda_matrix
# First, let's see if we can match them by position or need to subset beta_matrix

cat("\nIdentifying which 76 SVs were modeled...\n")

# Option 1: If the model dropped some SVs, we need to find which ones remain
# Let's check the mcmc.sp list length
cat("Number of species in mcmc.sp:", length(model_simple$mcmc.sp), "\n")

# The 76 modeled species should be the first 76 or a specific subset
# Let's create a proper subset of beta_matrix that matches lambda_matrix

# If beta_matrix has 405 rows but only 76 were modeled, 
# we need to subset beta_matrix to match
beta_matrix_modeled <- beta_matrix[1:76, ]  # Take first 76 rows
rownames(beta_matrix_modeled) <- rownames(lambda_matrix)  # Use lambda names

# Now recalculate top SVs from the MODELED set
top_n_svs <- 30
top_svs_idx_modeled <- order(rowSums(abs(beta_matrix_modeled)), decreasing = TRUE)[1:top_n_svs]

# Calculate residual correlation
lambda_cor <- lambda_matrix %*% t(lambda_matrix)
rownames(lambda_cor) <- rownames(lambda_matrix)
colnames(lambda_cor) <- rownames(lambda_matrix)

# Subset to top SVs
lambda_cor_subset <- lambda_cor[top_svs_idx_modeled, top_svs_idx_modeled]

pheatmap(lambda_cor_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Residual Correlation",
         fontsize_row = 8,
         fontsize_col = 8)

# For consistency, let's also redo Panels A and B with the modeled SVs
cat("\n--- Redoing all panels with the 76 modeled SVs ---\n")

# Panel A with modeled SVs
beta_subset_modeled <- beta_matrix_modeled[top_svs_idx_modeled, ]

pheatmap(beta_subset_modeled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-4, 4, length.out = 101),
         main = "A) Environmental Response (Modeled SVs)",
         fontsize_row = 8,
         fontsize_col = 10)

# Panel B with modeled SVs
sv_correlations_modeled <- cor(t(beta_matrix_modeled))
sv_cor_subset_modeled <- sv_correlations_modeled[top_svs_idx_modeled, top_svs_idx_modeled]

pheatmap(sv_cor_subset_modeled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "B) Covariance Matrix (Modeled SVs)",
         fontsize_row = 8,
         fontsize_col = 8)

# Panel C - need to subset pa_matrix to the 76 modeled SVs
pa_matrix_modeled <- pa_matrix[, 1:76]
colnames(pa_matrix_modeled) <- rownames(lambda_matrix)

sv_cooccurrence_modeled <- cor(pa_matrix_modeled)
sv_cooc_subset_modeled <- sv_cooccurrence_modeled[top_svs_idx_modeled, top_svs_idx_modeled]

pheatmap(sv_cooc_subset_modeled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "C) Co-occurrence in Data (Modeled SVs)",
         fontsize_row = 8,
         fontsize_col = 8)

# Step 21: Create visualizations of SV composition

# Add pond and family info for analysis
sv_richness <- data.frame(
  sample_id = rownames(pa_matrix),
  pond_id = env_data_scaled$pond_id,
  Family = env_data_scaled$Family,
  SV_richness = rowSums(pa_matrix)
)

# Summary statistics
cat("\nSV Richness Summary by Pond:\n")
pond_summary <- sv_richness %>% 
  group_by(pond_id) %>% 
  summarise(
    mean_richness = round(mean(SV_richness), 1),
    sd_richness = round(sd(SV_richness), 1),
    n_samples = n()
  )
print(pond_summary)

cat("\nSV Richness Summary by Family:\n")
family_summary <- sv_richness %>% 
  group_by(Family) %>% 
  summarise(
    mean_richness = round(mean(SV_richness), 1),
    sd_richness = round(sd(SV_richness), 1),
    n_samples = n()
  )
print(family_summary)

# Bar plot: SV richness by pond
plot_pond <- ggplot(sv_richness, aes(x = pond_id, y = SV_richness, fill = pond_id)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "SV Richness by Pond",
       x = "Pond ID",
       y = "Number of SVs Present") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_pond)

# Bar plot: SV richness by family
plot_family <- ggplot(sv_richness, aes(x = Family, y = SV_richness, fill = Family)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "SV Richness by Tadpole Family",
       x = "Family",
       y = "Number of SVs Present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_family)

# Combined plot: SV richness by family within each pond
plot_combined <- ggplot(sv_richness, aes(x = Family, y = SV_richness, fill = Family)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ pond_id, ncol = 5) +
  labs(title = "SV Richness by Family across Ponds",
       x = "Family",
       y = "Number of SVs Present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 9))

print(plot_combined)

