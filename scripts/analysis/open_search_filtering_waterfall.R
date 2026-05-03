# ============================================================
# Open Search Analysis Pipeline
# ============================================================
rm(list=ls())
# Load libraries
library(here)
library(tidyverse)

# ============================================================
# SET UP FOLDER STRUCTURE
# ============================================================

# Define output folders
csv_out <- here("data", "processed")
fig_out <- here("results", "figures", "open_search")

# Create folders if they don't already exist
dir.create(csv_out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================

# Load raw PTM candidates file
ptm_raw <- read_tsv(here("data", "processed", "PTM1_corrected.tsv"))

# Check it loaded correctly
nrow(ptm_raw)
colnames(ptm_raw)

# ============================================================
# STEP 1: FILTER TO KEY COMPARISONS
# ============================================================

key_comparisons <- c(
  "LPS / Unstim",
  "LPS+Nigericin / LPS",
  "LPS+Nigericin / Unstim",
  "P3C4 / Unstim",
  "P3C4+Nigericin / P3C4",
  "P3C4+Nigericin / Unstim"
)

ptm_filtered <- ptm_raw %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons)

# Check how many rows remain
nrow(ptm_filtered)

# ============================================================
# STEP 2: REMOVE CONTAMINANTS
# ============================================================

ptm_filtered <- ptm_filtered %>%
  filter(!grepl("CON__", UniProtIds))

# Check how many rows remain after contaminant removal
nrow(ptm_filtered)



# ============================================================
# STEP 3: FILTER MODIFICATIONS
# ============================================================
ptm_filtered %>%
  count(PTM.ModificationTitle) %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Define biological modifications to keep
biological_mods <- c(
  "Phospho (STY)",
  "Acetyl (K)",
  "GlyGly (K)",
  "Trimethyl (K)",
  "Sulfo (STY)",
  "Nitrosyl",
  "Crotonyl",
  "HexNAc (ST)",
  "Myristoyl",
  "Methyl (KR)",
  "Dimethyl (KR)",
  "Hex (K)"  # kept as reference artefact
)

# Filter to biological modifications only
ptm_bio <- ptm_filtered %>%
  filter(PTM.ModificationTitle %in% biological_mods)

# Check row count
nrow(ptm_bio)

# Check modification counts after filtering
ptm_bio %>%
  count(PTM.ModificationTitle) %>%
  arrange(desc(n))
# ============================================================
# STEP 4: PROTEIN ABUNDANCE CORRECTION
# ============================================================

# Load protein candidates file
protein_raw <- read_tsv(here("data_raw", "ptm_open_search_PG_candidates.tsv"))

# Flip comparisons in protein file
protein_flipped <- protein_raw %>%
  mutate(
    `Comparison (group1/group2)` = paste(
      `Condition Denominator`, "/", `Condition Numerator`
    ),
    `AVG Log2 Ratio` = -`AVG Log2 Ratio`
  )

# Filter protein file to key comparisons only
protein_filtered <- protein_flipped %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons)

# Check
nrow(protein_filtered)

# Join PTM file with protein file on UniProtIds and Comparison
# then subtract protein log2FC from PTM log2FC
ptm_corrected <- ptm_bio %>%
  left_join(
    protein_filtered %>%
      select(
        `Comparison (group1/group2)`,
        UniProtIds,
        protein_log2fc = `AVG Log2 Ratio`
      ),
    by = c("Comparison (group1/group2)", "UniProtIds")
  ) %>%
  mutate(
    corrected_log2fc = `AVG Log2 Ratio` - protein_log2fc
  )

# Check how many rows have a corrected value
ptm_corrected %>%
  summarise(
    total = n(),
    with_correction = sum(!is.na(corrected_log2fc)),
    without_correction = sum(is.na(corrected_log2fc))
  )
# ============================================================
# STEP 5: REMOVE UNCORRECTED ROWS AND DEFINE ROBUST SITES
# ============================================================

# Remove rows without protein correction
ptm_corrected <- ptm_corrected %>%
  filter(!is.na(corrected_log2fc))

# Check row count
nrow(ptm_corrected)

# Define robust sites
# Robust = p < 0.05 AND |corrected log2FC| >= 0.58
ptm_robust <- ptm_corrected %>%
  mutate(
    site_class = case_when(
      Pvalue < 0.05 & abs(corrected_log2fc) >= 0.58 ~ "Robust",
      Pvalue < 0.05 & abs(corrected_log2fc) < 0.58  ~ "Significant but small FC",
      Pvalue >= 0.05 & abs(corrected_log2fc) >= 0.58 ~ "Large FC but not significant",
      TRUE ~ "Not significant"
    )
  )

# Check how many robust sites
ptm_robust %>%
  count(site_class) %>%
  arrange(desc(n))

# Check robust sites by modification
ptm_robust %>%
  filter(site_class == "Robust") %>%
  count(PTM.ModificationTitle) %>%
  arrange(desc(n))
# ============================================================
# STEP 6: SAVE OUTPUTS
# ============================================================

# Save full corrected dataset
write_csv(
  ptm_corrected,
  file.path(csv_out, "ptm_open_search_corrected.csv")
)

# Save robust sites only
write_csv(
  ptm_robust %>% filter(site_class == "Robust"),
  file.path(csv_out, "ptm_open_search_robust.csv")
)

# Save all sites with classification
write_csv(
  ptm_robust,
  file.path(csv_out, "ptm_open_search_all_classified.csv")
)

message("Files saved successfully")
# ============================================================
# WATERFALL FIGURE
# ============================================================

# Define waterfall data manually using pipeline numbers
waterfall_data <- tibble(
  step = c(
    "Raw Spectronaut\nexport",
    "Contaminants\nremoved",
    "Key comparisons\nfilter",
    "Biological\nmodifications",
    "Abundance\ncorrected",
    "Robust sites"
  ),
  rows = c(828005, 809169, 331383, 58176, 58096, 1408)
) %>%
  mutate(
    step = factor(step, levels = step),
    percentage = round((rows / 828005) * 100, 1),
    label = paste0(format(rows, big.mark = ","), "\n(", percentage, "%)"),
    bar_colour = ifelse(step == "Robust sites", "Robust", "Other")
  )
# Plot
waterfall_plot <- ggplot(waterfall_data, aes(x = step, y = rows)) +
  geom_col(aes(fill = step == "Robust sites"), width = 0.6) +
  scale_fill_manual(
    values = c("FALSE" = "#2166AC", "TRUE" = "#D6191B"),
    guide = "none") +
  geom_text(
    aes(label = label),
    vjust = -0.4,
    size = 3.2,
    lineheight = 0.9
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "PTM candidates: filtering stages",
    x = NULL,
    y = "Number of rows retained"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.title.y = element_text(size = 10)
  )

# Display
print(waterfall_plot)

# Save
ggsave(
  filename = file.path(fig_out, "waterfall_filtering.png"),
  plot = waterfall_plot,
  width = 8,
  height = 5,
  dpi = 300
)

message("Waterfall figure saved")
# ============================================================
# PANEL B: Unique proteins per PTM class per comparison
# ============================================================

panel_b_data <- ptm_bio %>%
  group_by(`Comparison (group1/group2)`, PTM.ModificationTitle) %>%
  summarise(unique_proteins = n_distinct(UniProtIds), .groups = "drop")

panel_b <- ggplot(
  panel_b_data,
  aes(
    x = `Comparison (group1/group2)`,
    y = PTM.ModificationTitle,
    fill = unique_proteins
  )
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = unique_proteins), size = 3) +
  scale_fill_gradient(
    low = "#DEEBF7",
    high = "#2166AC",
    name = "Unique\nproteins"
  ) +
  labs(
    title = "Unique proteins per PTM class across key comparisons",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

print(panel_b)

ggsave(
  filename = file.path(fig_out, "panel_b_unique_proteins_heatmap.png"),
  plot = panel_b,
  width = 9,
  height = 6,
  dpi = 300
)
# ============================================================
# PANEL C: Mean absolute log2FC before vs after correction
# ============================================================

# Calculate mean absolute log2FC before and after correction
panel_c_data <- ptm_corrected %>%
  group_by(`Comparison (group1/group2)`, PTM.ModificationTitle) %>%
  summarise(
    before = mean(abs(`AVG Log2 Ratio`), na.rm = TRUE),
    after = mean(abs(corrected_log2fc), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(before, after),
    names_to = "correction",
    values_to = "mean_abs_log2fc"
  ) %>%
  mutate(
    correction = factor(
      correction,
      levels = c("before", "after"),
      labels = c("Before correction", "After correction")
    )
  )

# Plot
panel_c <- ggplot(
  panel_c_data,
  aes(
    x = `Comparison (group1/group2)`,
    y = PTM.ModificationTitle,
    fill = mean_abs_log2fc
  )
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean_abs_log2fc, 2)), size = 2.8) +
  scale_fill_gradient(
    low = "#FFF5EB",
    high = "#D6191B",
    name = "Mean |log2FC|"
  ) +
  facet_wrap(~correction) +
  labs(
    title = "Mean absolute log2FC: before vs after protein abundance correction",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    strip.background = element_rect(fill = "grey85", colour = "grey50"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

print(panel_c)

ggsave(
  filename = file.path(fig_out, "panel_c_fc_correction_heatmap.png"),
  plot = panel_c,
  width = 12,
  height = 6,
  dpi = 300
)
# ============================================================
# PANEL D: Robust sites per PTM class per comparison
# ============================================================

panel_d_data <- ptm_robust %>%
  filter(site_class == "Robust") %>%
  group_by(`Comparison (group1/group2)`, PTM.ModificationTitle) %>%
  summarise(robust_sites = n(), .groups = "drop") %>%
  # Fill in zeros for modification/comparison combos with no robust sites
  complete(
    `Comparison (group1/group2)`,
    PTM.ModificationTitle,
    fill = list(robust_sites = 0)
  )

panel_d <- ggplot(
  panel_d_data,
  aes(
    x = `Comparison (group1/group2)`,
    y = PTM.ModificationTitle,
    fill = robust_sites
  )
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = robust_sites), size = 3) +
  scale_fill_gradient(
    low = "#EFF3FF",
    high = "#2166AC",
    name = "Robust\nsites"
  ) +
  labs(
    title = "Robust PTM sites per modification across key comparisons",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

print(panel_d)

ggsave(
  filename = file.path(fig_out, "panel_d_robust_sites_heatmap.png"),
  plot = panel_d,
  width = 9,
  height = 6,
  dpi = 300
)
