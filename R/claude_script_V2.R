# =============================================================================
#  DCNF Joint Species Distribution Model — Top 50 SVs
#  Replicates the Python analysis in R / RStudio
#
#  Required packages (install once):
#    install.packages(c("tidyverse", "vegan", "pheatmap", "RColorBrewer",
#                       "ggplot2", "patchwork", "scales", "ggrepel"))
#
#  Set your working directory to the folder containing the CSV files, or
#  update the file paths below.
# =============================================================================


# ── 0. Libraries ──────────────────────────────────────────────────────────────

library(tidyverse)      # data wrangling + ggplot2
library(vegan)          # ecology multivariate tools (PCA/RDA, Bray-Curtis, etc.)
library(pheatmap)       # clustered heatmap
library(RColorBrewer)   # colour palettes
library(patchwork)      # combine ggplots
library(scales)         # percent / comma formatting
library(ggrepel)        # non-overlapping labels on ordination


# ── 1. Load data ──────────────────────────────────────────────────────────────

# Adjust paths if your files are elsewhere
meta <- read_csv("data/raw/DCNF_metadata.csv", show_col_types = FALSE)
otu  <- read_csv("data/raw/otu_table.csv",     show_col_types = FALSE)
tax  <- read_csv("data/raw/tax_table.csv",     show_col_types = FALSE)


# ── 2. Identify DCNF samples ──────────────────────────────────────────────────
#
#  The OTU table uses abbreviated pond codes (no dashes):
#    08  181  191  291  292  294  361  871  SFA2  SFA4

pond_codes <- c("08", "181", "191", "291", "292", "294",
                "361", "871", "SFA2", "SFA4")

# Function: extract pond code from sample name
get_pond <- function(s) {
  for (p in pond_codes) {
    if (grepl(p, s, fixed = TRUE)) return(p)
  }
  return(NA_character_)
}

# Function: extract sample type from sample name prefix
get_type <- function(s) {
  prefixes <- c("WSUB", "WWAT", "RC", "SI", "AT", "GC",
                "GT", "BF", "AC", "NV", "RS")
  for (t in prefixes) {
    if (startsWith(s, t)) return(t)
  }
  return("other")
}

# Build sample metadata from OTU sample names
sample_info <- tibble(SampleName = otu$SampleName) %>%
  mutate(
    Pond = map_chr(SampleName, get_pond),
    Type = map_chr(SampleName, get_type)
  ) %>%
  filter(!is.na(Pond))          # keep only DCNF samples

dcnf_samples <- sample_info$SampleName

cat("DCNF samples found:", length(dcnf_samples), "\n")
cat("Ponds:", paste(unique(sample_info$Pond), collapse = ", "), "\n")
cat("Sample types:", paste(unique(sample_info$Type), collapse = ", "), "\n")


# ── 3. Subset OTU table to DCNF samples ───────────────────────────────────────

otu_dcnf <- otu %>%
  filter(SampleName %in% dcnf_samples) %>%
  column_to_rownames("SampleName")    # rows = samples, cols = SVs


# ── 4. Top 50 most abundant SVs ───────────────────────────────────────────────

sv_totals <- colSums(otu_dcnf)
top50_svs <- names(sort(sv_totals, decreasing = TRUE))[1:50]

otu_top50 <- otu_dcnf[, top50_svs]

cat("\nTop 5 SVs by total count:\n")
print(head(sort(sv_totals[top50_svs], decreasing = TRUE), 5))


# ── 5. Taxonomy labels ────────────────────────────────────────────────────────

tax_top50 <- tax %>%
  filter(SV %in% top50_svs) %>%
  select(SV, Phylum, Class, Order, Family, Genus) %>%
  slice(match(top50_svs, SV))          # keep top50 order

# Short label: first non-NA level from Genus → Family → Order → Class → Phylum
make_label <- function(sv) {
  row <- tax_top50 %>% filter(SV == sv)
  for (col in c("Genus", "Family", "Order", "Class", "Phylum")) {
    val <- row[[col]]
    if (!is.na(val) && nchar(trimws(val)) > 0 && val != "NA") {
      return(paste0(val, "\n(", sv, ")"))
    }
  }
  return(sv)
}

sv_labels <- setNames(map_chr(top50_svs, make_label), top50_svs)


# ── 6. Transformations ────────────────────────────────────────────────────────

# 6a. Relative abundance (%)
otu_ra <- otu_top50 / rowSums(otu_top50) * 100

# 6b. Log(count + 1) — for heatmap
otu_log <- log1p(otu_top50)

# 6c. Hellinger — for ordination
otu_hell <- decostand(otu_ra / 100, method = "hellinger")


# ── 7. Colour palettes ────────────────────────────────────────────────────────

pond_order <- sort(unique(sample_info$Pond))
type_order  <- sort(unique(sample_info$Type))

pond_cols <- setNames(brewer.pal(max(3, length(pond_order)), "Set2")[seq_along(pond_order)],
                      pond_order)
type_cols  <- setNames(brewer.pal(max(3, length(type_order)), "Set1")[seq_along(type_order)],
                       type_order)


# =============================================================================
#  FIGURE 1 — Clustered Heatmap (pheatmap)
# =============================================================================

# Annotation data frame (rows = samples)
ann_row <- sample_info %>%
  column_to_rownames("SampleName") %>%
  select(Pond, Type)

ann_colors <- list(
  Pond = pond_cols,
  Type = type_cols
)

# Rename columns with taxonomy labels
mat_heat          <- as.matrix(otu_log)
colnames(mat_heat) <- sv_labels[colnames(mat_heat)]

png("fig1_heatmap.png", width = 2600, height = 1800, res = 120)
pheatmap(
  mat_heat,
  annotation_row    = ann_row,
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("white", "#FEB24C", "#BD0026"))(100),
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  fontsize_row      = 5,
  fontsize_col      = 6,
  angle_col         = 90,
  main              = "DCNF — Top 50 SVs (Ward clustering, log counts)",
  border_color      = NA
)
dev.off()
cat("Fig 1 saved\n")


# =============================================================================
#  FIGURE 2 — Mean Relative Abundance by Pond (stacked bar)
# =============================================================================

# Join relative abundances with pond info, compute means per pond per SV
pond_mean_long <- otu_ra %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  pivot_longer(cols = all_of(top50_svs),
               names_to = "SV", values_to = "RA") %>%
  group_by(Pond, SV) %>%
  summarise(Mean_RA = mean(RA), .groups = "drop") %>%
  mutate(Label = sv_labels[SV])

# 20-colour palette for 50 SVs (cycling tab20)
sv_fill_cols <- setNames(
  rep(c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Accent")),
      length.out = 50),
  top50_svs
)

fig2 <- ggplot(pond_mean_long, aes(x = Pond, y = Mean_RA, fill = SV)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = sv_fill_cols,
                    labels = sv_labels,
                    guide  = guide_legend(ncol = 1, keyheight = unit(0.35, "cm"))) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Mean Relative Abundance of Top 50 SVs by Pond (DCNF)",
    x = "Pond", y = "Mean Relative Abundance (%)", fill = "SV (Taxonomy)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        axis.text.x  = element_text(size = 10))

ggsave("fig2_pond_bar.png", fig2, width = 22, height = 9, dpi = 150)
cat("Fig 2 saved\n")


# =============================================================================
#  FIGURE 3 — Per-Sample Stacked Bar (sorted by Pond → Type)
# =============================================================================

sample_order <- sample_info %>%
  arrange(Pond, Type) %>%
  pull(SampleName)

per_sample_long <- otu_ra %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  pivot_longer(cols = all_of(top50_svs),
               names_to = "SV", values_to = "RA") %>%
  mutate(
    SampleName = factor(SampleName, levels = sample_order),
    SampleLabel = paste(Pond, Type, sep = "/")
  )

fig3 <- ggplot(per_sample_long,
               aes(x = SampleName, y = RA, fill = SV)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_manual(values = sv_fill_cols, labels = sv_labels,
                    guide = guide_legend(ncol = 1, keyheight = unit(0.3, "cm"))) +
  labs(
    title = "Per-Sample Relative Abundance of Top 50 SVs — DCNF",
    x = "Sample (sorted by Pond / Type)",
    y = "Relative Abundance (%)",
    fill = "SV (Taxonomy)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, size = 5),
    legend.text  = element_text(size = 4.5),
    legend.title = element_text(size = 6)
  )

# Add vertical dashed lines at pond boundaries
pond_breaks <- sample_info %>%
  arrange(Pond, Type) %>%
  mutate(idx = row_number()) %>%
  group_by(Pond) %>%
  summarise(start = min(idx), .groups = "drop") %>%
  filter(start > 1) %>%
  pull(start)

fig3 <- fig3 +
  geom_vline(xintercept = pond_breaks - 0.5,
             linetype = "dashed", colour = "black", linewidth = 0.7)

ggsave("fig3_sample_bar.png", fig3, width = 32, height = 9, dpi = 150)
cat("Fig 3 saved\n")


# =============================================================================
#  FIGURE 4 — PCA Ordination (Hellinger-transformed)
# =============================================================================

pca_res    <- rda(otu_hell)            # vegan's rda() = PCA when no constraints
pca_scores <- as.data.frame(scores(pca_res, display = "sites", scaling = 1))
var_exp    <- round(eigenvals(pca_res)[1:2] / sum(eigenvals(pca_res)) * 100, 1)

pca_df <- pca_scores %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  rename(PC1 = PC1, PC2 = PC2)   # vegan uses PC1/PC2 naming

make_pca_plot <- function(df, colour_var, pal, title_suffix) {
  ggplot(df, aes(x = PC1, y = PC2,
                 colour = .data[[colour_var]],
                 shape  = .data[[colour_var]])) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 8, 4, 3, 7), 4)[seq_along(pal)]) +
    labs(
      title  = paste("PCA — coloured by", title_suffix),
      x      = paste0("PC1 (", var_exp[1], "%)"),
      y      = paste0("PC2 (", var_exp[2], "%)"),
      colour = title_suffix,
      shape  = title_suffix
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
}

fig4a <- make_pca_plot(pca_df, "Pond", pond_cols, "Pond")
fig4b <- make_pca_plot(pca_df, "Type", type_cols, "Sample Type")

fig4 <- (fig4a | fig4b) +
  plot_annotation(
    title = "DCNF — PCA Ordination (Hellinger-transformed relative abundance)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("fig4_pca.png", fig4, width = 16, height = 7, dpi = 150)
cat("Fig 4 saved\n")


# =============================================================================
#  FIGURE 5 — Boxplots of Top 12 SVs across Ponds
# =============================================================================

top12 <- top50_svs[1:12]

box_long <- otu_ra[, top12] %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  pivot_longer(cols = all_of(top12),
               names_to = "SV", values_to = "RA") %>%
  mutate(
    Label = sv_labels[SV],
    Label = factor(Label, levels = sv_labels[top12]),  # preserve rank order
    Pond  = factor(Pond, levels = pond_order)
  )

fig5 <- ggplot(box_long, aes(x = Pond, y = RA, fill = Pond)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ Label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pond_cols) +
  labs(
    title = "Relative Abundance Distribution of Top 12 SVs by Pond (DCNF)",
    x     = "Pond",
    y     = "Relative Abundance (%)",
    fill  = "Pond"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
    strip.text    = element_text(size = 6.5),
    legend.position = "bottom"
  )

ggsave("fig5_boxplots.png", fig5, width = 20, height = 14, dpi = 150)
cat("Fig 5 saved\n")


# =============================================================================
#  BONUS — Bray-Curtis distance heatmap between samples
# =============================================================================

bc_dist <- vegdist(otu_ra / 100, method = "bray")
bc_mat  <- as.matrix(bc_dist)

png("fig6_braycurtis.png", width = 1800, height = 1600, res = 120)
pheatmap(
  bc_mat,
  annotation_row    = ann_row,
  annotation_col    = ann_row,
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
  clustering_method = "ward.D2",
  fontsize_row      = 5,
  fontsize_col      = 5,
  main              = "Bray-Curtis Dissimilarity Between DCNF Samples (Top 50 SVs)",
  border_color      = NA
)
dev.off()
cat("Fig 6 (Bray-Curtis heatmap) saved\n")


# =============================================================================
#  BONUS — PERMANOVA: do ponds explain community composition?
# =============================================================================

# Match sample order between distance matrix and metadata
meta_ordered <- sample_info[match(rownames(bc_mat), sample_info$SampleName), ]

set.seed(42)
permanova_pond <- adonis2(
  bc_dist ~ Pond,
  data   = meta_ordered,
  permutations = 999
)
cat("\n── PERMANOVA: Pond effect ──\n")
print(permanova_pond)

permanova_type <- adonis2(
  bc_dist ~ Type,
  data   = meta_ordered,
  permutations = 999
)
cat("\n── PERMANOVA: Sample Type effect ──\n")
print(permanova_type)

permanova_both <- adonis2(
  bc_dist ~ Pond + Type,
  data   = meta_ordered,
  permutations = 999
)
cat("\n── PERMANOVA: Pond + Type ──\n")
print(permanova_both)


# =============================================================================
#  BONUS — NMDS ordination (alternative to PCA)
# =============================================================================

set.seed(42)
nmds_res <- metaMDS(otu_ra / 100, distance = "bray", k = 2,
                    trymax = 100, trace = FALSE)
cat("\nNMDS stress:", round(nmds_res$stress, 4), "\n")

nmds_df <- as.data.frame(scores(nmds_res, display = "sites")) %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName")

fig_nmds_pond <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2,
                                     colour = Pond, shape = Pond)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = pond_cols) +
  scale_shape_manual(values = rep(c(16,17,15,18,8,4,3,7), 2)[seq_along(pond_cols)]) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("Stress = ", round(nmds_res$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 3, colour = "grey40") +
  labs(title = "NMDS — coloured by Pond",
       x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 11)

fig_nmds_type <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2,
                                     colour = Type, shape = Type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = type_cols) +
  scale_shape_manual(values = rep(c(16,17,15,18,8,4,3,7), 2)[seq_along(type_cols)]) +
  labs(title = "NMDS — coloured by Sample Type",
       x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 11)

fig_nmds <- (fig_nmds_pond | fig_nmds_type) +
  plot_annotation(
    title = "DCNF — NMDS Ordination (Bray-Curtis, Top 50 SVs)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("fig7_nmds.png", fig_nmds, width = 16, height = 7, dpi = 150)
cat("Fig 7 (NMDS) saved\n")


# =============================================================================
#  Summary table — top 50 SVs with taxonomy and total counts
# =============================================================================

summary_tbl <- tax_top50 %>%
  mutate(Total_reads = sv_totals[SV],
         Label       = sv_labels[SV]) %>%
  arrange(desc(Total_reads)) %>%
  select(SV, Label, Phylum, Class, Order, Family, Genus, Total_reads)

write_csv(summary_tbl, "DCNF_top50_SVs_summary.csv")
cat("\nSummary table saved → DCNF_top50_SVs_summary.csv\n")
cat("\nAll done!\n")
