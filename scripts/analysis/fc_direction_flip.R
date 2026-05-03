rm (list = ls())
ptm_biological <- fread("data/processed/PTM1_biological.csv")
protein_da     <- fread("data_raw/ptm_open_search_PG_candidates.tsv")
# =============================================================
# Step 2 addition: Subset to biological comparisons, 
# flip log2FC direction where needed
# Applies to BOTH ptm_biological and protein_da
# =============================================================

# --- Define comparisons to EXCLUDE (ambiguous) ---
exclude_comparisons <- c(
  "LPS / Nigericin",
  "P3C4 / Nigericin",
  "P3C4 / LPS"
)

# --- Define comparisons to FLIP ---
# Currently less stimulated / more stimulated
# Need to become more stimulated / less stimulated
flip_comparisons <- c(
  "Unstim / LPS",
  "Unstim / P3C4",
  "Unstim / Nigericin",
  "Unstim / LPS+Nigericin",
  "Unstim / P3C4+Nigericin",
  "LPS / LPS+Nigericin",
  "P3C4 / P3C4+Nigericin",
  "P3C4 / LPS+Nigericin"
)

# --- Define comparisons to KEEP AS IS ---
keep_comparisons <- c(
  "P3C4+Nigericin / LPS",
  "P3C4+Nigericin / Nigericin",
  "P3C4+Nigericin / LPS+Nigericin",
  "LPS+Nigericin / Nigericin"
)

# Verify all 15 comparisons are accounted for
all_comparisons <- unique(ptm_biological$`Comparison (group1/group2)`)
accounted_for   <- c(exclude_comparisons, 
                     flip_comparisons, 
                     keep_comparisons)

cat("Total comparisons:", length(all_comparisons), "\n")
cat("Accounted for:", length(accounted_for), "\n")
cat("Any missing?\n")
print(setdiff(all_comparisons, accounted_for))

# =============================================================
# Apply to PTM biological dataset
# =============================================================

ptm_corrected_direction <- ptm_biological %>%
  # Remove ambiguous comparisons
  filter(!`Comparison (group1/group2)` %in% exclude_comparisons) %>%
  # Flip log2FC and swap numerator/denominator labels
  mutate(
    needs_flip = `Comparison (group1/group2)` %in% flip_comparisons,
    
    # Flip ratio values
    `AVG Log2 Ratio` = ifelse(needs_flip,
                              -`AVG Log2 Ratio`,
                              `AVG Log2 Ratio`),
    `Absolute AVG Log2 Ratio` = abs(`AVG Log2 Ratio`),
    `% Change` = ifelse(needs_flip,
                        -`% Change`,
                        `% Change`),
    
    # Swap numerator and denominator labels
    new_numerator   = ifelse(needs_flip,
                             `Condition Denominator`,
                             `Condition Numerator`),
    new_denominator = ifelse(needs_flip,
                             `Condition Numerator`,
                             `Condition Denominator`),
    
    # Update comparison label to reflect new direction
    `Comparison (group1/group2)` = ifelse(
      needs_flip,
      paste(new_numerator, "/", new_denominator),
      `Comparison (group1/group2)`
    ),
    `Condition Numerator`   = new_numerator,
    `Condition Denominator` = new_denominator
  ) %>%
  select(-needs_flip, -new_numerator, -new_denominator)

# =============================================================
# Apply same logic to protein DA
# =============================================================

protein_da_corrected <- protein_da %>%
  filter(!`Comparison (group1/group2)` %in% exclude_comparisons) %>%
  mutate(
    needs_flip = `Comparison (group1/group2)` %in% flip_comparisons,
    
    `AVG Log2 Ratio` = ifelse(needs_flip,
                              -`AVG Log2 Ratio`,
                              `AVG Log2 Ratio`),
    `Absolute AVG Log2 Ratio` = abs(`AVG Log2 Ratio`),
    `% Change` = ifelse(needs_flip,
                        -`% Change`,
                        `% Change`),
    
    new_numerator   = ifelse(needs_flip,
                             `Condition Denominator`,
                             `Condition Numerator`),
    new_denominator = ifelse(needs_flip,
                             `Condition Numerator`,
                             `Condition Denominator`),
    
    `Comparison (group1/group2)` = ifelse(
      needs_flip,
      paste(new_numerator, "/", new_denominator),
      `Comparison (group1/group2)`
    ),
    `Condition Numerator`   = new_numerator,
    `Condition Denominator` = new_denominator
  ) %>%
  select(-needs_flip, -new_numerator, -new_denominator)

# =============================================================
# Verify the flip worked correctly
# =============================================================

cat("\nComparisons remaining in PTM dataset:\n")
print(unique(ptm_corrected_direction$`Comparison (group1/group2)`))

cat("\nComparisons remaining in protein DA dataset:\n")
print(unique(protein_da_corrected$`Comparison (group1/group2)`))

# Spot check - pick one flipped comparison and confirm
# log2FC direction has changed
cat("\nSpot check LPS+Nigericin / LPS (was LPS / LPS+Nigericin):\n")
cat("Mean AVG Log2 Ratio after flip:\n")
print(mean(ptm_corrected_direction$`AVG Log2 Ratio`[
  ptm_corrected_direction$`Comparison (group1/group2)` == 
    "LPS+Nigericin / LPS"], na.rm = TRUE))

# Should be opposite sign to original
cat("Mean AVG Log2 Ratio in original (LPS / LPS+Nigericin):\n")
print(mean(ptm_biological$`AVG Log2 Ratio`[
  ptm_biological$`Comparison (group1/group2)` == 
    "LPS / LPS+Nigericin"], na.rm = TRUE))

# =============================================================
# Export
# =============================================================

write.csv(ptm_corrected_direction,
          "data/processed/PTM1_biological_corrected.csv",
          row.names = FALSE)

write.csv(protein_da_corrected,
          "data/processed/protein_da_corrected.csv",
          row.names = FALSE)

message("Direction correction complete.")
message("PTM rows retained: ", nrow(ptm_corrected_direction))
message("Protein DA rows retained: ", nrow(protein_da_corrected))