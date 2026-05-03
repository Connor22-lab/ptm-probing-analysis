# =============================================================
# Step 3: Protein abundance correction
# Run this as a standalone script from scratch
# Inputs:
#   data/processed/PTM1_biological_corrected.csv
#   data/processed/protein_da_corrected.csv
# =============================================================
rm(list = ls())
library(dplyr)
library(data.table)

# --- 1. Load files ---
ptm <- fread("data/processed/PTM1_biological_corrected.csv")
protein_da <- fread("data/processed/protein_da_corrected.csv")

# --- 2. Clean and deduplicate protein_da ---
protein_da_clean <- protein_da %>%
  # Remove contaminants
  filter(!grepl("^Cont", Group)) %>%
  # Rename AVG Log2 Ratio to distinguish from PTM ratio
  rename(Protein_AVG_Log2Ratio = `AVG Log2 Ratio`) %>%
  # Remove duplicates - keep row with most ratios
  group_by(Genes, `Comparison (group1/group2)`) %>%
  slice_max(order_by = `# of Ratios`,
            n = 1,
            with_ties = FALSE) %>%
  ungroup()

cat("Protein DA rows after cleaning:", 
    nrow(protein_da_clean), "\n")

# --- 3. Merge ---
ptm_corrected <- ptm %>%
  left_join(
    protein_da_clean %>%
      select(Genes,
             `Comparison (group1/group2)`,
             Protein_AVG_Log2Ratio),
    by = c("Genes", "Comparison (group1/group2)"),
    relationship = "many-to-one"
  )

cat("Rows before merge:", nrow(ptm), "\n")
cat("Rows after merge:", nrow(ptm_corrected), "\n")
cat("Row count unchanged:", 
    nrow(ptm) == nrow(ptm_corrected), "\n")
cat("NAs in Protein_AVG_Log2Ratio:", 
    sum(is.na(ptm_corrected$Protein_AVG_Log2Ratio)), "\n")

# Check NAs
cat("NA percentage:", 
    round(1637 / 88360 * 100, 1), "%\n")

# Which genes are unmatched
unmatched_genes <- ptm_corrected %>%
  filter(is.na(Protein_AVG_Log2Ratio)) %>%
  distinct(Genes) %>%
  pull(Genes)

cat("Unique unmatched genes:", 
    length(unmatched_genes), "\n")

# Are any POI in unmatched
poi_list <- c("NLRP3", "PYCARD", "CASP1", "GSDMD",
              "IL1B", "IL18", "HSP90AA1", "HSP90AB1",
              "HSPA8", "PPP2CA")

cat("Any POI unmatched:\n")
print(intersect(unmatched_genes, poi_list))
# --- 4. Calculate corrected log2FC ---
ptm_corrected <- ptm_corrected %>%
  mutate(
    PTM_log2FC           = `AVG Log2 Ratio`,
    Corrected_log2FC     = PTM_log2FC - Protein_AVG_Log2Ratio,
    Abs_corrected_log2FC = abs(Corrected_log2FC),
    Correction_class = case_when(
      is.na(Protein_AVG_Log2Ratio) ~ "No protein data",
      abs(Protein_AVG_Log2Ratio) < 0.3 ~ "Protein stable",
      sign(Corrected_log2FC) != sign(PTM_log2FC) ~
        "Direction reversed by correction",
      Abs_corrected_log2FC < abs(PTM_log2FC) * 0.5 ~
        "Substantially reduced by correction",
      TRUE ~ "Robust after correction"
    )
  )

cat("\nCorrection classification summary:\n")
print(table(ptm_corrected$Correction_class))

# --- 5. PTM class selection table with corrected values ---
key_comparisons <- c(
  "LPS / Unstim",
  "P3C4 / Unstim",
  "LPS+Nigericin / LPS",
  "P3C4+Nigericin / P3C4",
  "LPS+Nigericin / Unstim",
  "P3C4+Nigericin / Unstim"
)

# --- Fixed PTM class selection table ---
selection_corrected <- ptm_corrected %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons,
         !is.na(Corrected_log2FC)) %>%
  group_by(PTM.ModificationTitle,
           `Comparison (group1/group2)`) %>%
  summarise(
    Mean_abs_PTM_log2FC  = mean(abs(PTM_log2FC), 
                                na.rm = TRUE),
    Mean_abs_corr_log2FC = mean(Abs_corrected_log2FC,
                                na.rm = TRUE),
    N_unique_proteins    = n_distinct(Genes),
    N_robust             = sum(Correction_class == 
                                 "Robust after correction"),
    N_total              = n(),
    .groups = "drop"
  ) %>%
  group_by(PTM.ModificationTitle) %>%
  summarise(
    Max_mean_abs_PTM_log2FC  = round(max(Mean_abs_PTM_log2FC), 3),
    Max_mean_abs_corr_log2FC = round(max(Mean_abs_corr_log2FC), 3),
    Mean_proteins            = round(mean(N_unique_proteins), 0),
    Pct_robust               = round(sum(N_robust) / 
                                       sum(N_total) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_mean_abs_corr_log2FC))

print(selection_corrected)
write.csv(selection_corrected,
          "results/tables/step3_ptm_selection_corrected.csv",
          row.names = FALSE)

# --- 6. POI corrected summary ---
poi_list <- c("NLRP3", "PYCARD", "CASP1", "GSDMD",
              "IL1B", "IL18", "HSP90AA1", "HSP90AB1",
              "HSPA8", "PPP2CA")

poi_corrected <- ptm_corrected %>%
  filter(Genes %in% poi_list,
         `Comparison (group1/group2)` %in% key_comparisons) %>%
  group_by(Genes, PTM.ModificationTitle,
           `Comparison (group1/group2)`) %>%
  summarise(
    N_sites              = n_distinct(Group),
    PTM_log2FC           = round(mean(PTM_log2FC,
                                      na.rm = TRUE), 3),
    Protein_log2FC       = round(mean(Protein_AVG_Log2Ratio,
                                      na.rm = TRUE), 3),
    Corrected_log2FC     = round(mean(Corrected_log2FC,
                                      na.rm = TRUE), 3),
    Abs_corrected_log2FC = round(mean(Abs_corrected_log2FC,
                                      na.rm = TRUE), 3),
    Pvalue               = round(min(Pvalue, na.rm = TRUE), 4),
    Correction_class     = first(Correction_class),
    .groups = "drop"
  ) %>%
  arrange(Genes, desc(Abs_corrected_log2FC))

print(poi_corrected)
write.csv(poi_corrected,
          "results/tables/step3_poi_corrected.csv",
          row.names = FALSE)

# --- 7. Export ---
write.csv(ptm_corrected,
          "data/processed/PTM1_abundance_corrected.csv",
          row.names = FALSE)

message("Step 3 complete.")