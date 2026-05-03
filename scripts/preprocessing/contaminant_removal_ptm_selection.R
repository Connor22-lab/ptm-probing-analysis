# =============================================================
# Step 2 final: Contaminant removal + biological PTM selection
# Input:  data_raw/PTM1.tsv (already localisation-filtered 
#         at source by Spectronaut at 0.75 threshold)
# =============================================================
rm(list = ls())
library(dplyr)
library(data.table)
library(ggplot2)

ptm_candidates <- fread("data_raw/PTM1.tsv")

# --- Stage 1: Remove contaminants ---
n_raw <- nrow(ptm_candidates)

ptm_nocont <- ptm_candidates %>%
  mutate(UniProt_extracted = sub("_.*", "", Group)) %>%
  filter(!grepl("^Cont", UniProt_extracted))

n_nocont <- nrow(ptm_nocont)
cat("Contaminants removed:", n_raw - n_nocont, "\n")

# --- Stage 2: Retain biological PTMs only ---
biological_ptms <- c(
  "Phospho (STY)",
  "GlyGly (K)",
  "Acetyl (K)",
  "HexNAc (ST)",
  "Trimethyl (K)",
  "Crotonyl",
  "Myristoyl",
  "Hex (K)",
  "Sulfo (STY)",
  "Nitrosyl",
  "Oxidation (W)"
)

ptm_biological <- ptm_nocont %>%
  filter(PTM.ModificationTitle %in% biological_ptms)

n_biological <- nrow(ptm_biological)
cat("Rows after biological PTM filter:", n_biological, "\n")
cat("Unique proteins:", n_distinct(ptm_biological$Genes), "\n")
cat("Unique sites:", n_distinct(ptm_biological$Group), "\n")

# --- Summary by modification type ---
mod_summary <- ptm_biological %>%
  group_by(PTM.ModificationTitle) %>%
  summarise(
    N_rows            = n(),
    N_unique_proteins = n_distinct(Genes),
    N_unique_sites    = n_distinct(Group),
    N_comparisons     = n_distinct(`Comparison (group1/group2)`),
    .groups = "drop"
  ) %>%
  arrange(desc(N_unique_proteins))

print(mod_summary)
write.csv(mod_summary,
          "results/tables/step2_biological_ptm_summary.csv",
          row.names = FALSE)

# --- Filtering waterfall ---
filter_summary <- data.frame(
  Stage = c(
    "Raw (post-Spectronaut export)",
    "Contaminants removed",
    "Biological PTMs only"
  ),
  N_rows = c(n_raw, n_nocont, n_biological)
) %>%
  mutate(
    N_removed    = c(0,
                     n_raw    - n_nocont,
                     n_nocont - n_biological),
    Pct_retained = round(N_rows / n_raw * 100, 1),
    Stage        = factor(Stage, levels = Stage)
  )

print(filter_summary)
write.csv(filter_summary,
          "results/tables/step2_filtering_summary_final.csv",
          row.names = FALSE)

p_waterfall <- ggplot(filter_summary,
                      aes(x = Stage, y = N_rows)) +
  geom_bar(stat = "identity", fill = "#2C7BB6",
           colour = "black", width = 0.6) +
  geom_text(aes(label = paste0(N_rows, "\n(",
                               Pct_retained, "%)")),
            vjust = -0.3, size = 3) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "PTM candidates: filtering stages",
       x = NULL,
       y = "Number of rows retained") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("results/figures/step2_filter_waterfall.png", p_waterfall,
       width = 8, height = 5, dpi = 300)

# --- Export ---
write.csv(ptm_biological,
          "data/processed/PTM1_biological.csv",
          row.names = FALSE)

message("Step 2 complete.")
