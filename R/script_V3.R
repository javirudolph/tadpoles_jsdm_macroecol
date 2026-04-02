# =============================================================================
#  DCNF Joint Species Distribution Model — Top 50 SVs, classified by ORDER
#
#  SVs are labelled and grouped at the Order level only (no finer taxonomy).
#  Figures produced:
#    fig1_heatmap.png        — Clustered heatmap with Order colour bar on SVs
#    fig2_pond_bar.png       — Mean relative abundance by Pond, stacked by Order
#    fig3_sample_bar.png     — Per-sample stacked bar, stacked by Order
#    fig4_pca.png            — PCA ordination coloured by Pond and Sample Type
#    fig5_order_boxplots.png — Boxplots of every Order's RA across Ponds
#    fig6_braycurtis.png     — Bray-Curtis dissimilarity heatmap
#    fig7_nmds.png           — NMDS ordination
#    DCNF_top50_SVs_summary.csv
#
#  Required packages (install once):
#    install.packages(c("tidyverse", "vegan", "pheatmap",
#                       "RColorBrewer", "patchwork", "scales"))
#
#  Place all CSV files in your working directory (Session → Set Working
#  Directory → To Source File Location is the easiest approach in RStudio).
# =============================================================================


# ── 0. Libraries ──────────────────────────────────────────────────────────────

library(tidyverse)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(scales)


# ── 1. Load data ──────────────────────────────────────────────────────────────

meta <- read_csv("DCNF_metadata.csv", show_col_types = FALSE)
otu  <- read_csv("otu_table.csv",     show_col_types = FALSE)
tax  <- read_csv("tax_table.csv",     show_col_types = FALSE)


# ── 2. Identify DCNF samples in the OTU table ─────────────────────────────────
#
#  The OTU table uses abbreviated, dash-free pond codes embedded in sample
#  names: 08  181  191  291  292  294  361  871  SFA2  SFA4

pond_codes <- c("08", "181", "191", "291", "292", "294",
                "361", "871", "SFA2", "SFA4")

get_pond <- function(s) {
  for (p in pond_codes) {
    if (grepl(p, s, fixed = TRUE)) return(p)
  }
  return(NA_character_)
}

get_type <- function(s) {
  for (t in c("WSUB", "WWAT", "RC", "SI", "AT", "GC",
              "GT", "BF", "AC", "NV", "RS")) {
    if (startsWith(s, t)) return(t)
  }
  return("other")
}

sample_info <- tibble(SampleName = otu$SampleName) %>%
  mutate(
    Pond = map_chr(SampleName, get_pond),
    Type = map_chr(SampleName, get_type)
  ) %>%
  filter(!is.na(Pond))

dcnf_samples <- sample_info$SampleName

cat("DCNF samples found :", length(dcnf_samples), "\n")
cat("Ponds              :", paste(sort(unique(sample_info$Pond)), collapse = ", "), "\n")
cat("Sample types       :", paste(sort(unique(sample_info$Type)), collapse = ", "), "\n")


# ── 3. Subset OTU table ───────────────────────────────────────────────────────

otu_dcnf <- otu %>%
  filter(SampleName %in% dcnf_samples) %>%
  column_to_rownames("SampleName")      # rows = samples, cols = SVs


# ── 4. Top 50 most abundant SVs ───────────────────────────────────────────────

sv_totals <- colSums(otu_dcnf)
top50_svs <- names(sort(sv_totals, decreasing = TRUE))[1:50]
otu_top50 <- otu_dcnf[, top50_svs]

cat("\nTop 10 SVs by total count:\n")
print(head(sort(sv_totals[top50_svs], decreasing = TRUE), 10))


# ── 5. Taxonomy at Order level only ───────────────────────────────────────────

tax_top50 <- tax %>%
  filter(SV %in% top50_svs) %>%
  select(SV, Phylum, Class, Order) %>%    # stop at Order
  slice(match(top50_svs, SV))             # preserve rank order

# Verify every SV has an Order assigned
stopifnot(all(!is.na(tax_top50$Order)))
cat("\nAll top-50 SVs have an Order assigned ✓\n")
cat("Unique Orders:", tax_top50$Order %>% unique() %>% length(), "\n")
cat(paste(sort(unique(tax_top50$Order)), collapse = "\n  "), "\n")

# SV label: "Order\n(SVxx)"
sv_labels <- setNames(
  paste0(tax_top50$Order, "\n(", tax_top50$SV, ")"),
  tax_top50$SV
)

# Named vector: SV → Order
sv_order_map <- setNames(tax_top50$Order, tax_top50$SV)


# ── 6. Transformations ────────────────────────────────────────────────────────

# Relative abundance (%)
otu_ra <- otu_top50 / rowSums(otu_top50) * 100

# Log(count + 1) for heatmap
otu_log <- log1p(otu_top50)

# Hellinger for ordination
otu_hell <- decostand(otu_ra / 100, method = "hellinger")


# ── 7. Order-level aggregated relative abundance ──────────────────────────────
#
#  Collapse the 50 SVs into 18 Orders by summing relative abundances.
#  This is the key step that differs from the Genus-level script.

# Rename SV columns to their Order, then sum columns sharing the same Order
otu_ra_ord <- otu_ra
colnames(otu_ra_ord) <- sv_order_map[colnames(otu_ra_ord)]

# sapply trick: for each unique order, rowSums of matching columns
order_list <- sort(unique(sv_order_map))

otu_ra_by_order <- as.data.frame(
  sapply(order_list, function(ord) {
    cols <- which(colnames(otu_ra_ord) == ord)
    if (length(cols) == 1) otu_ra_ord[, cols]
    else rowSums(otu_ra_ord[, cols])
  })
)
# rows = samples, cols = orders; sums to ~100% per sample
# (small deviations due to only using top 50 SVs)


# ── 8. Colour palettes ────────────────────────────────────────────────────────

# 18 qualitative colours for Orders — combined from tab20 equivalents in R
order_pal_raw <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
  "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
  "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
  "#c49c94","#f7b6d2","#c7c7c7"
)
order_cols <- setNames(order_pal_raw[seq_along(order_list)], order_list)

# Sort orders by total abundance across all DCNF samples (most abundant first)
order_totals_vec <- colSums(otu_ra_by_order)
order_sorted <- names(sort(order_totals_vec, decreasing = TRUE))

# Pond and type palettes
pond_order_vec <- sort(unique(sample_info$Pond))
type_order_vec  <- sort(unique(sample_info$Type))

pond_cols <- setNames(
  brewer.pal(max(3, length(pond_order_vec)), "Set2")[seq_along(pond_order_vec)],
  pond_order_vec
)
type_cols <- setNames(
  brewer.pal(max(3, length(type_order_vec)), "Set1")[seq_along(type_order_vec)],
  type_order_vec
)


# =============================================================================
#  FIGURE 1 — Clustered heatmap with Order colour bar on SV columns
# =============================================================================

# Row annotation (samples)
ann_row <- sample_info %>%
  column_to_rownames("SampleName") %>%
  select(Pond, Type)

# Column annotation (SVs → their Order)
ann_col <- data.frame(
  Order = sv_order_map[top50_svs],
  row.names = sv_labels[top50_svs]    # label is "Order\n(SVxx)"
)

ann_colors <- list(
  Pond  = pond_cols,
  Type  = type_cols,
  Order = order_cols
)

mat_heat           <- as.matrix(otu_log)
colnames(mat_heat) <- sv_labels[colnames(mat_heat)]

png("fig1_heatmap.png", width = 3000, height = 1900, res = 120)
pheatmap(
  mat_heat,
  annotation_row    = ann_row,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("white", "#FEB24C", "#BD0026"))(100),
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  fontsize_row      = 5,
  fontsize_col      = 6,
  angle_col         = 90,
  main              = "DCNF — Top 50 SVs classified by Order (Ward clustering, log counts)",
  border_color      = NA
)
dev.off()
cat("Fig 1 saved → fig1_heatmap.png\n")


# =============================================================================
#  FIGURE 2 — Mean relative abundance by Pond, stacked by Order
# =============================================================================

pond_mean_ord <- otu_ra_by_order %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  pivot_longer(cols = all_of(order_list),
               names_to = "Order", values_to = "RA") %>%
  group_by(Pond, Order) %>%
  summarise(Mean_RA = mean(RA), .groups = "drop") %>%
  mutate(Order = factor(Order, levels = order_sorted))   # most abundant first

fig2 <- ggplot(pond_mean_ord,
               aes(x = Pond, y = Mean_RA, fill = Order)) +
  geom_bar(stat = "identity", width = 0.78) +
  scale_fill_manual(values = order_cols,
                    breaks = order_sorted,          # legend sorted by abundance
                    guide = guide_legend(ncol = 1,
                                        keyheight = unit(0.45, "cm"))) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Mean Relative Abundance of Top 50 SVs by Pond — grouped by Order (DCNF)",
    x = "Pond", y = "Mean Relative Abundance (%)", fill = "Order"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    axis.text.x  = element_text(size = 10)
  )

ggsave("fig2_pond_bar.png", fig2, width = 16, height = 9, dpi = 150)
cat("Fig 2 saved → fig2_pond_bar.png\n")


# =============================================================================
#  FIGURE 3 — Per-sample stacked bar (sorted Pond → Type)
# =============================================================================

sample_order_vec <- sample_info %>%
  arrange(Pond, Type) %>%
  pull(SampleName)

per_sample_ord <- otu_ra_by_order %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  pivot_longer(cols = all_of(order_list),
               names_to = "Order", values_to = "RA") %>%
  mutate(
    SampleName = factor(SampleName, levels = sample_order_vec),
    Order      = factor(Order, levels = order_sorted)
  )

# Pond boundary positions (for vertical separator lines)
pond_breaks <- sample_info %>%
  arrange(Pond, Type) %>%
  mutate(idx = row_number()) %>%
  group_by(Pond) %>%
  summarise(start = min(idx), .groups = "drop") %>%
  filter(start > 1) %>%
  pull(start)

fig3 <- ggplot(per_sample_ord,
               aes(x = SampleName, y = RA, fill = Order)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_vline(xintercept = pond_breaks - 0.5,
             linetype = "dashed", colour = "black", linewidth = 0.8) +
  scale_fill_manual(values = order_cols,
                    breaks = order_sorted,
                    guide = guide_legend(ncol = 1,
                                        keyheight = unit(0.45, "cm"))) +
  labs(
    title = "Per-Sample Relative Abundance of Top 50 SVs — grouped by Order (DCNF)",
    x     = "Sample (sorted by Pond → Type)",
    y     = "Relative Abundance (%)",
    fill  = "Order"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, size = 5),
    legend.text  = element_text(size = 7.5),
    legend.title = element_text(size = 8.5)
  )

ggsave("fig3_sample_bar.png", fig3, width = 32, height = 9, dpi = 150)
cat("Fig 3 saved → fig3_sample_bar.png\n")


# =============================================================================
#  FIGURE 4 — PCA Ordination (Hellinger-transformed)
# =============================================================================

pca_res  <- rda(otu_hell)
pca_sc   <- as.data.frame(scores(pca_res, display = "sites", scaling = 1))
var_exp  <- round(eigenvals(pca_res)[1:2] / sum(eigenvals(pca_res)) * 100, 1)

pca_df <- pca_sc %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName")

make_pca_plot <- function(df, colour_var, pal, title_sfx) {
  ggplot(df, aes(x = PC1, y = PC2,
                 colour = .data[[colour_var]],
                 shape  = .data[[colour_var]])) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(
      values = rep(c(16,17,15,18,8,4,3,7), 4)[seq_along(pal)]
    ) +
    labs(
      title  = paste("PCA — coloured by", title_sfx),
      x      = paste0("PC1 (", var_exp[1], "%)"),
      y      = paste0("PC2 (", var_exp[2], "%)"),
      colour = title_sfx, shape = title_sfx
    ) +
    theme_bw(base_size = 11)
}

fig4 <- (make_pca_plot(pca_df, "Pond", pond_cols, "Pond") |
         make_pca_plot(pca_df, "Type", type_cols, "Sample Type")) +
  plot_annotation(
    title = "DCNF — PCA Ordination (Hellinger-transformed, Top 50 SVs)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("fig4_pca.png", fig4, width = 17, height = 7, dpi = 150)
cat("Fig 4 saved → fig4_pca.png\n")


# =============================================================================
#  FIGURE 5 — Boxplots: every Order's RA distribution across Ponds
# =============================================================================

box_long_ord <- otu_ra_by_order %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName") %>%
  pivot_longer(cols = all_of(order_list),
               names_to = "Order", values_to = "RA") %>%
  mutate(
    Order = factor(Order, levels = order_sorted),   # most abundant first
    Pond  = factor(Pond,  levels = pond_order_vec)
  )

fig5 <- ggplot(box_long_ord, aes(x = Pond, y = RA, fill = Pond)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
  facet_wrap(~ Order, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pond_cols) +
  labs(
    title = "Relative Abundance (%) of All 18 Orders Across Ponds — DCNF Top 50 SVs",
    x = "Pond", y = "Relative Abundance (%)", fill = "Pond"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    strip.text      = element_text(size = 7.5, face = "bold"),
    legend.position = "bottom"
  )

ggsave("fig5_order_boxplots.png", fig5, width = 22, height = 16, dpi = 150)
cat("Fig 5 saved → fig5_order_boxplots.png\n")


# =============================================================================
#  FIGURE 6 — Bray-Curtis dissimilarity heatmap
# =============================================================================

bc_dist <- vegdist(otu_ra / 100, method = "bray")
bc_mat  <- as.matrix(bc_dist)

png("fig6_braycurtis.png", width = 1900, height = 1700, res = 120)
pheatmap(
  bc_mat,
  annotation_row    = ann_row,
  annotation_col    = ann_row,
  annotation_colors = ann_colors[c("Pond", "Type")],
  color             = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
  clustering_method = "ward.D2",
  fontsize_row      = 5,
  fontsize_col      = 5,
  main              = "Bray-Curtis Dissimilarity — DCNF Samples (Top 50 SVs)",
  border_color      = NA
)
dev.off()
cat("Fig 6 saved → fig6_braycurtis.png\n")


# =============================================================================
#  FIGURE 7 — NMDS ordination (Bray-Curtis)
# =============================================================================

set.seed(42)
nmds_res <- metaMDS(otu_ra / 100, distance = "bray", k = 2,
                    trymax = 100, trace = FALSE)
cat("\nNMDS stress:", round(nmds_res$stress, 4), "\n")

nmds_df <- as.data.frame(scores(nmds_res, display = "sites")) %>%
  rownames_to_column("SampleName") %>%
  left_join(sample_info, by = "SampleName")

stress_label <- paste0("Stress = ", round(nmds_res$stress, 3))

fig7_pond <- ggplot(nmds_df,
                    aes(x = NMDS1, y = NMDS2, colour = Pond, shape = Pond)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = pond_cols) +
  scale_shape_manual(
    values = rep(c(16,17,15,18,8,4,3,7), 2)[seq_along(pond_cols)]
  ) +
  annotate("text", x = Inf, y = -Inf, label = stress_label,
           hjust = 1.1, vjust = -0.5, size = 3, colour = "grey40") +
  labs(title = "NMDS — coloured by Pond",
       x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 11)

fig7_type <- ggplot(nmds_df,
                    aes(x = NMDS1, y = NMDS2, colour = Type, shape = Type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = type_cols) +
  scale_shape_manual(
    values = rep(c(16,17,15,18,8,4,3,7), 2)[seq_along(type_cols)]
  ) +
  labs(title = "NMDS — coloured by Sample Type",
       x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 11)

fig7 <- (fig7_pond | fig7_type) +
  plot_annotation(
    title = "DCNF — NMDS Ordination (Bray-Curtis, Top 50 SVs)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("fig7_nmds.png", fig7, width = 16, height = 7, dpi = 150)
cat("Fig 7 saved → fig7_nmds.png\n")


# =============================================================================
#  PERMANOVA — do ponds and sample types explain community structure?
# =============================================================================

meta_ordered <- sample_info[match(rownames(bc_mat), sample_info$SampleName), ]

set.seed(42)
cat("\n── PERMANOVA: Pond ──────────────────────────────────\n")
print(adonis2(bc_dist ~ Pond, data = meta_ordered, permutations = 999))

cat("\n── PERMANOVA: Sample Type ───────────────────────────\n")
print(adonis2(bc_dist ~ Type, data = meta_ordered, permutations = 999))

cat("\n── PERMANOVA: Pond + Sample Type ────────────────────\n")
print(adonis2(bc_dist ~ Pond + Type, data = meta_ordered, permutations = 999))


# =============================================================================
#  Summary table
# =============================================================================

summary_tbl <- tax_top50 %>%
  mutate(
    Total_reads = sv_totals[SV],
    SV_label    = sv_labels[SV]
  ) %>%
  arrange(desc(Total_reads)) %>%
  select(SV, SV_label, Phylum, Class, Order, Total_reads)

write_csv(summary_tbl, "DCNF_top50_SVs_summary.csv")
cat("\nSummary table saved → DCNF_top50_SVs_summary.csv\n")
cat("\n✓ All done!\n")
