# =============================================================
# Standalone FIGURE SCRIPT (3 figures)
# Inputs:
#   1) tables/ptm_class_rank.csv        (from your standalone prioritisation script)
#   2) tables/protein_hotspots.csv      (from your standalone prioritisation script)
#   3) data/processed/PTM1_final_corrected.csv  (site-level corrected dataset)
#
# Outputs (PNG):
#   figures/Fig1_PTM_class_prioritisation_bubble.png
#   figures/Fig2_Protein_hotspots_bar.png
#   figures/Fig3_PTMxComparison_robust_heatmap.png
# =============================================================

rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)

# -------------------------
# Files + thresholds
# -------------------------
ptm_rank_file   <- "results/tables/ptm_class_rank.csv"
hotspots_file   <- "results/tables/protein_hotspots.csv"
ptm_final_file  <- "data/processed/PTM1_final_corrected.csv"

# Must match how you defined "robust" in your pipeline
fc_thresh <- 0.58
p_thresh  <- 0.05

# Your key comparisons (set order on heatmap)
key_comparisons <- c(
  "LPS / Unstim",
  "P3C4 / Unstim",
  "LPS+Nigericin / LPS",
  "P3C4+Nigericin / P3C4",
  "LPS+Nigericin / Unstim",
  "P3C4+Nigericin / Unstim"
)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load tables
# -------------------------
ptm_rank <- fread(ptm_rank_file) %>% as_tibble()
hotspots <- fread(hotspots_file) %>% as_tibble()

# -------------------------
# FIGURE 1: PTM class prioritisation (bubble)
#   x = Mean_abs_corr
#   size = Frac_robust_sites
#   colour = Direction_consistency
#   y ordered by Composite_score
# -------------------------
ptm_rank_plot <- ptm_rank %>%
  mutate(
    PTM.ModificationTitle = factor(
      PTM.ModificationTitle,
      levels = PTM.ModificationTitle[order(Composite_score, decreasing = TRUE)]
    ),
    Frac_robust_pct = 100 * Frac_robust_sites
  )

p1 <- ggplot(ptm_rank_plot,
             aes(x = Mean_abs_corr,
                 y = PTM.ModificationTitle,
                 size = Frac_robust_pct,
                 colour = Direction_consistency)) +
  geom_point(alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(
    title = "PTM class prioritisation after protein abundance correction",
    x = "Mean |Corrected log2FC|",
    y = NULL,
    size = "Robust sites (%)",
    colour = "Direction\nconsistency"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

ggsave("results/figures/Fig1_PTM_class_prioritisation_bubble.png",
       p1, width = 9.5, height = 5.8, dpi = 300)

# -------------------------
# FIGURE 2: Protein hotspots (bar)
#   x = N_robust_sites
#   y = Genes (top N)
#   fill = PTM_classes
# -------------------------
topN <- 25

hotspots_plot <- hotspots %>%
  slice_head(n = topN) %>%
  mutate(
    Genes = factor(Genes, levels = rev(Genes)),  # keep ranking order top->bottom
    PTM_classes = as.factor(PTM_classes)
  )

p2 <- ggplot(hotspots_plot,
             aes(x = N_robust_sites, y = Genes, fill = PTM_classes)) +
  geom_col(width = 0.75) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Top ", topN, " proteins with robust PTM-specific regulation"),
    x = "Number of robust PTM sites (p < 0.05 and |corr log2FC| ≥ 0.58)",
    y = NULL,
    fill = "PTM\nclasses"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 9)
  )

ggsave("results/figures/Fig2_Protein_hotspots_bar.png",
       p2, width = 10.5, height = 7.2, dpi = 300)

# -------------------------
# FIGURE 3: Condition specificity heatmap
#   PTM class × key comparisons
#   Fill = number of robust sites
# -------------------------
# -------------------------
# FIGURE 3: Condition specificity heatmap
#   PTM class × key comparisons
#   Fill = FRACTION of robust sites
# -------------------------

ptm_final <- fread(ptm_final_file) %>% as_tibble()

# Define robust site at site-level
ptm_final <- ptm_final %>%
  mutate(
    Robust_site = !is.na(Corrected_log2FC) &
      !is.na(Pvalue) &
      (Pvalue < p_thresh) &
      (abs(Corrected_log2FC) >= fc_thresh)
  ) %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons)

# Use top PTM classes (by composite score)
top_ptmN <- 10
top_ptms <- ptm_rank %>%
  arrange(desc(Composite_score)) %>%
  slice_head(n = top_ptmN) %>%
  pull(PTM.ModificationTitle)

heatmap_df <- ptm_final %>%
  filter(PTM.ModificationTitle %in% top_ptms) %>%
  group_by(PTM.ModificationTitle, `Comparison (group1/group2)`) %>%
  summarise(
    N_total_sites  = n(),
    N_robust_sites = sum(Robust_site, na.rm = TRUE),
    Frac_robust = ifelse(N_total_sites > 0,
                         N_robust_sites / N_total_sites,
                         0),
    .groups = "drop"
  ) %>%
  mutate(
    `Comparison (group1/group2)` = factor(`Comparison (group1/group2)`,
                                          levels = key_comparisons),
    PTM.ModificationTitle = factor(
      PTM.ModificationTitle,
      levels = top_ptms
    )
  )

p3 <- ggplot(heatmap_df,
             aes(x = `Comparison (group1/group2)`,
                 y = PTM.ModificationTitle,
                 fill = Frac_robust)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", Frac_robust)),
            size = 3) +
  scale_fill_gradient(low = "white", high = "#D73027") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    title = paste0("Fraction of robust PTM sites by condition (Top ", top_ptmN, " PTM classes)"),
    x = NULL,
    y = NULL,
    fill = "Fraction\nrobust"
  )

ggsave("results/figures/Fig3_PTMxComparison_fraction_robust_heatmap.png",
       p3, width = 12, height = 6.5, dpi = 300)

cat("\nSaved figures to /figures:\n")
cat("  - Fig1_PTM_class_prioritisation_bubble.png\n")
cat("  - Fig2_Protein_hotspots_bar.png\n")
cat("  - Fig3_PTMxComparison_robust_heatmap.png\n\n")
