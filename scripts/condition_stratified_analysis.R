# =============================================================
# Standalone Step 3: PTM abundance correction + data-driven prioritisation
# =============================================================
rm(list = ls())

library(data.table)
library(dplyr)

# -------------------------
# User parameters
# -------------------------
key_comparisons <- c(
  "LPS / Unstim",
  "P3C4 / Unstim",
  "LPS+Nigericin / LPS",
  "P3C4+Nigericin / P3C4",
  "LPS+Nigericin / Unstim",
  "P3C4+Nigericin / Unstim"
)

fc_thresh <- 0.58   # effect size threshold (log2 scale)
p_thresh  <- 0.05   # p-value threshold

# Optional: keep your current POI list, but we won't rely on it
poi_list <- c("NLRP3", "PYCARD", "CASP1", "GSDMD",
              "IL1B", "IL18", "HSP90AA1", "HSP90AB1",
              "HSPA8", "PPP2CA")

# -------------------------
# I/O
# -------------------------
ptm_file     <- "data_processed/PTM1_biological_corrected.csv"
protein_file <- "data_processed/protein_da_corrected.csv"

dir.create("tables", showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load
# -------------------------
ptm_raw <- fread(ptm_file)
protein_da <- fread(protein_file)

# -------------------------
# Clean + deduplicate protein DA (one row per Genes × Comparison)
# -------------------------
protein_da_clean <- protein_da %>%
  filter(!grepl("^Cont", Group)) %>%
  rename(Protein_AVG_Log2Ratio = `AVG Log2 Ratio`) %>%
  group_by(Genes, `Comparison (group1/group2)`) %>%
  slice_max(order_by = `# of Ratios`, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Genes, `Comparison (group1/group2)`, Protein_AVG_Log2Ratio)

# -------------------------
# Merge + compute corrected values
# -------------------------
ptm <- ptm_raw %>%
  left_join(protein_da_clean,
            by = c("Genes", "Comparison (group1/group2)")) %>%
  mutate(
    PTM_log2FC = `AVG Log2 Ratio`,
    Corrected_log2FC = PTM_log2FC - Protein_AVG_Log2Ratio,
    Abs_corrected_log2FC = abs(Corrected_log2FC),
    
    # Site-level robust flag (simple + defensible)
    Robust_site = !is.na(Corrected_log2FC) &
      !is.na(Pvalue) &
      (Pvalue < p_thresh) &
      (Abs_corrected_log2FC >= fc_thresh)
  )

# Filter to comparisons of interest
ptm_key <- ptm %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons)

# -------------------------
# PTM class ranking
# -------------------------
ptm_class_rank <- ptm_key %>%
  filter(!is.na(Corrected_log2FC)) %>%
  group_by(PTM.ModificationTitle) %>%
  summarise(
    N_sites = n(),
    N_proteins = n_distinct(Genes),
    Mean_abs_corr = mean(Abs_corrected_log2FC, na.rm = TRUE),
    Frac_robust_sites = mean(Robust_site, na.rm = TRUE),
    
    # Optional: direction consistency (1 = mostly one direction; 0 = mixed)
    Direction_consistency = abs(mean(sign(Corrected_log2FC), na.rm = TRUE)),
    
    # A simple composite score (you can report it or ignore it)
    Composite_score = Mean_abs_corr * Frac_robust_sites * Direction_consistency,
    .groups = "drop"
  ) %>%
  arrange(desc(Composite_score))

write.csv(ptm_class_rank, "tables/ptm_class_rank.csv", row.names = FALSE)

# -------------------------
# Protein hotspot ranking (data-driven)
# -------------------------
protein_hotspots <- ptm_key %>%
  filter(Robust_site) %>%
  group_by(Genes) %>%
  summarise(
    N_robust_sites = n(),
    PTM_classes = n_distinct(PTM.ModificationTitle),
    Comparisons = n_distinct(`Comparison (group1/group2)`),
    Max_abs_corr = max(Abs_corrected_log2FC, na.rm = TRUE),
    Mean_abs_corr = mean(Abs_corrected_log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(PTM_classes), desc(Comparisons), desc(N_robust_sites), desc(Max_abs_corr))

write.csv(protein_hotspots, "tables/protein_hotspots.csv", row.names = FALSE)

# -------------------------
# Optional: POI summary (computed correctly; no "first()")
# -------------------------
poi_summary <- ptm_key %>%
  filter(Genes %in% poi_list, !is.na(Corrected_log2FC)) %>%
  group_by(Genes, PTM.ModificationTitle, `Comparison (group1/group2)`) %>%
  summarise(
    N_sites = n_distinct(Group),
    
    # Important: these are different summaries (keep both)
    Mean_corr = mean(Corrected_log2FC, na.rm = TRUE),
    Mean_abs_corr = mean(abs(Corrected_log2FC), na.rm = TRUE),
    
    N_robust_sites = sum(Robust_site, na.rm = TRUE),
    Min_p = min(Pvalue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(N_robust_sites), desc(Mean_abs_corr))

write.csv(poi_summary, "tables/poi_composite_summary.csv", row.names = FALSE)

# -------------------------
# Print a tiny preview so you know it's working
# -------------------------
cat("\nDone.\n")
cat("\nTop PTM classes:\n")
print(head(ptm_class_rank, 10))

cat("\nTop protein hotspots:\n")
print(head(protein_hotspots, 15))

cat("\nPOI summary (top rows):\n")
print(head(poi_summary, 15))

