rm(list = ls())
setwd("C:/Users/conno/Desktop/Final Year School/Final Year Project")
library(tidyverse)

# Load PTM data
ptm_data <- read_tsv("data_raw/PTM1.tsv")

# Define explicit mapping of what should be flipped
flip_mapping <- tribble(
  ~Original_Comparison,              ~Corrected_Comparison,
  "Unstim / LPS",                    "LPS / Unstim",
  "Unstim / P3C4",                   "P3C4 / Unstim",
  "Unstim / Nigericin",              "Nigericin / Unstim",
  "Unstim / LPS+Nigericin",          "LPS+Nigericin / Unstim",
  "Unstim / P3C4+Nigericin",         "P3C4+Nigericin / Unstim",
  "LPS / LPS+Nigericin",             "LPS+Nigericin / LPS",
  "P3C4 / P3C4+Nigericin",           "P3C4+Nigericin / P3C4"
)

# Join with mapping and correct
ptm_corrected <- ptm_data %>%
  left_join(flip_mapping, by = c("Comparison (group1/group2)" = "Original_Comparison")) %>%
  mutate(
    # Determine if this comparison needs flipping
    needs_flip = !is.na(Corrected_Comparison),
    
    # Use corrected comparison if available, otherwise keep original
    `Comparison (group1/group2)` = coalesce(Corrected_Comparison, `Comparison (group1/group2)`),
    
    # Extract new numerator and denominator from corrected comparison
    New_Numerator = str_split_fixed(`Comparison (group1/group2)`, " / ", 2)[,1],
    New_Denominator = str_split_fixed(`Comparison (group1/group2)`, " / ", 2)[,2],
    
    # Update condition columns (adjust column names as needed)
    `Condition (Numerator)` = New_Numerator,
    `Condition (Denominator)` = New_Denominator,
    
    # Flip log2FC if needed
    `AVG Log2 Ratio` = if_else(needs_flip, -`AVG Log2 Ratio`, `AVG Log2 Ratio`),
    `Absolute AVG Log2 Ratio` = abs(`AVG Log2 Ratio`)
  ) %>%
  # Clean up temporary columns
  select(-Corrected_Comparison, -needs_flip, -New_Numerator, -New_Denominator)

# Save
write_tsv(ptm_corrected, "data_processed/PTM1_corrected.tsv")

# Verify
print("Unique corrected comparisons:")
print(unique(ptm_corrected$`Comparison (group1/group2)`))
