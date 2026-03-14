# =============================================================
# Step 3: Condition-stratified analysis
# =============================================================
rm(list = ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# --- Load corrected datasets ---
ptm <- fread("data_processed/PTM1_biological_corrected.csv")
protein_da <- fread("data_processed/protein_da_corrected.csv")

# Confirm comparisons look correct
cat("Comparisons in corrected PTM dataset:\n")
print(unique(ptm$`Comparison (group1/group2)`))

# =============================================================
# 3A: PTM class overview across comparisons
# How many unique proteins and sites per PTM class 
# per comparison
# =============================================================

ptm_class_overview <- ptm %>%
  group_by(`Comparison (group1/group2)`, 
           PTM.ModificationTitle) %>%
  summarise(
    N_unique_proteins = n_distinct(Genes),
    N_unique_sites    = n_distinct(Group),
    N_rows            = n(),
    Mean_abs_log2FC   = round(mean(`Absolute AVG Log2 Ratio`,
                                   na.rm = TRUE), 3),
    Median_abs_log2FC = round(median(`Absolute AVG Log2 Ratio`,
                                     na.rm = TRUE), 3),
    .groups = "drop"
  )

write.csv(ptm_class_overview,
          "tables/step3_ptm_class_overview.csv",
          row.names = FALSE)

# --- Heatmap: N unique proteins per PTM class per comparison ---
p_heatmap_proteins <- ggplot(ptm_class_overview,
                             aes(x = `Comparison (group1/group2)`,
                                 y = PTM.ModificationTitle,
                                 fill = N_unique_proteins)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = N_unique_proteins), 
            size = 2.5) +
  scale_fill_gradientn(
    colours = c("white", "#2C7BB6", "#00008B"),
    name = "N proteins") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9)
  ) +
  labs(title = "PTM class coverage: unique proteins per comparison",
       x = NULL, y = NULL)

ggsave("figures/step3_heatmap_proteins_per_ptmclass.png",
       p_heatmap_proteins,
       width = 10, height = 7, dpi = 300)

# --- Heatmap: mean absolute log2FC per PTM class 
#     per comparison ---
p_heatmap_fc <- ggplot(ptm_class_overview,
                       aes(x = `Comparison (group1/group2)`,
                           y = PTM.ModificationTitle,
                           fill = Mean_abs_log2FC)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = round(Mean_abs_log2FC, 2)),
            size = 2.5) +
  scale_fill_gradientn(
    colours = c("white", "#FEB24C", "#F03B20"),
    name = "Mean |log2FC|") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9)
  ) +
  labs(title = "PTM class: mean absolute log2FC per comparison",
       x = NULL, y = NULL)

ggsave("figures/step3_heatmap_fc_per_ptmclass.png",
       p_heatmap_fc,
       width = 10, height = 7, dpi = 300)

# =============================================================
# 3B: Focus on key biological comparisons
# Priming: LPS / Unstim, P3C4 / Unstim
# Activation: LPS+Nigericin / LPS, P3C4+Nigericin / P3C4
# Fully activated vs baseline: LPS+Nigericin / Unstim,
#                               P3C4+Nigericin / Unstim
# =============================================================

key_comparisons <- c(
  "LPS / Unstim",
  "P3C4 / Unstim",
  "LPS+Nigericin / LPS",
  "P3C4+Nigericin / P3C4",
  "LPS+Nigericin / Unstim",
  "P3C4+Nigericin / Unstim"
)

ptm_key <- ptm %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons)

cat("\nRows in key comparisons:", nrow(ptm_key), "\n")

# --- Bar chart: N unique proteins per PTM class 
#     in key comparisons ---
ptm_key_summary <- ptm_key %>%
  group_by(PTM.ModificationTitle,
           `Comparison (group1/group2)`) %>%
  summarise(
    N_unique_proteins = n_distinct(Genes),
    N_unique_sites    = n_distinct(Group),
    Mean_abs_log2FC   = mean(`Absolute AVG Log2 Ratio`,
                             na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(`Comparison (group1/group2)` = factor(
    `Comparison (group1/group2)`, levels = key_comparisons))

p_bar_proteins <- ggplot(ptm_key_summary,
                         aes(x = reorder(PTM.ModificationTitle, N_unique_proteins),
                             y = N_unique_proteins,
                             fill = `Comparison (group1/group2)`)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 11) +
  labs(title = "Unique proteins per PTM class: key comparisons",
       x = NULL,
       y = "Number of unique proteins",
       fill = "Comparison")

ggsave("figures/step3_bar_proteins_key_comparisons.png",
       p_bar_proteins,
       width = 10, height = 7, dpi = 300)

# =============================================================
# 3C: Inflammasome protein focus
# Define POI list - expand beyond canonical NLRP3 components
# =============================================================

# Canonical NLRP3 inflammasome components
nlrp3_core <- c("NLRP3", "PYCARD", "CASP1", 
                "GSDMD", "IL1B", "IL18")

# Regulatory proteins - ubiquitin system
ubiquitin_reg <- c("BRCC3", "FBXL2", "FBXO3", 
                   "TRIM31", "MARCH7")

# Regulatory proteins - kinases/phosphatases
kinase_reg <- c("PKA", "JNK1", "MAPK8", "BTK", 
                "TBK1", "IKBKE", "PPP2CA")

# HSP chaperones - known NLRP3 regulators
chaperones <- c("HSP90AA1", "HSP90AB1", "HSPA1A", 
                "HSPA8")

# Other inflammasome sensors
other_sensors <- c("AIM2", "NLRC4", "PYRIN", "MEFV",
                   "NLRP1", "NLRP6")

# Combine all POI
poi_list <- c(nlrp3_core, ubiquitin_reg, 
              kinase_reg, chaperones, other_sensors)

# Filter PTM dataset to POI
ptm_poi <- ptm_key %>%
  filter(Genes %in% poi_list)

cat("\nPOI detected in dataset:\n")
print(unique(ptm_poi$Genes))
cat("Total POI rows:", nrow(ptm_poi), "\n")

# PTM coverage per POI protein
poi_ptm_coverage <- ptm_poi %>%
  group_by(Genes, PTM.ModificationTitle) %>%
  summarise(
    N_sites           = n_distinct(Group),
    N_comparisons     = n_distinct(`Comparison (group1/group2)`),
    Mean_abs_log2FC   = round(mean(`Absolute AVG Log2 Ratio`,
                                   na.rm = TRUE), 3),
    Max_abs_log2FC    = round(max(`Absolute AVG Log2 Ratio`,
                                  na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(Genes, desc(N_sites))

print(poi_ptm_coverage)
write.csv(poi_ptm_coverage,
          "tables/step3_poi_ptm_coverage.csv",
          row.names = FALSE)

# =============================================================
# 3D: Volcano plots for key comparisons
# Coloured by PTM class, POI labelled
# =============================================================

library(ggrepel)

# Function to generate one volcano per comparison
volcano_by_comparison <- function(data, comparison_name, 
                                  poi, fc_thresh = 0.58,
                                  p_thresh = 0.05) {
  
  df <- data %>%
    filter(`Comparison (group1/group2)` == comparison_name) %>%
    mutate(
      log2FC   = `AVG Log2 Ratio`,
      neg_logP = -log10(Pvalue),
      is_poi   = Genes %in% poi,
      label    = ifelse(is_poi, Genes, NA)
    )
  
  p <- ggplot(df, aes(x = log2FC, y = neg_logP,
                      colour = PTM.ModificationTitle)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(p_thresh),
               linetype = "dashed", colour = "grey40") +
    geom_point(data = filter(df, is_poi),
               aes(x = log2FC, y = neg_logP),
               colour = "black", size = 2.5, shape = 21,
               fill = NA, stroke = 1) +
    geom_label_repel(aes(label = label),
                     na.rm = TRUE,
                     size = 2.5,
                     colour = "black",
                     max.overlaps = 20,
                     box.padding = 0.3) +
    scale_colour_brewer(palette = "Set3",
                        name = "PTM class") +
    theme_bw(base_size = 11) +
    labs(
      title = comparison_name,
      x     = "AVG Log2 Ratio",
      y     = "-log10(p-value)"
    )
  
  return(p)
}

# Generate volcano for each key comparison
for (comp in key_comparisons) {
  p_vol <- volcano_by_comparison(
    data            = ptm_key,
    comparison_name = comp,
    poi             = poi_list
  )
  
  # Clean filename
  fname <- gsub("[/ ]", "_", comp)
  fname <- gsub("\\+", "plus", fname)
  
  ggsave(paste0("figures/step3_volcano_", fname, ".png"),
         p_vol,
         width = 9, height = 6, dpi = 300)
  
  message("Saved volcano: ", comp)
}

message("Step 3 complete.")

# Detailed look at PTMs on detected POI
poi_detail <- ptm_poi %>%
  group_by(Genes, PTM.ModificationTitle,
           `Comparison (group1/group2)`) %>%
  summarise(
    N_sites         = n_distinct(Group),
    Mean_log2FC     = round(mean(`AVG Log2 Ratio`, 
                                 na.rm = TRUE), 3),
    Mean_abs_log2FC = round(mean(`Absolute AVG Log2 Ratio`,
                                 na.rm = TRUE), 3),
    Min_Pvalue      = round(min(Pvalue, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(Genes, PTM.ModificationTitle)

print(poi_detail)
write.csv(poi_detail,
          "tables/step3_poi_ptm_detail.csv",
          row.names = FALSE)

# Also get a simple presence/absence summary
# Which PTM types does each POI carry
poi_presence <- ptm_poi %>%
  group_by(Genes, PTM.ModificationTitle) %>%
  summarise(
    N_sites       = n_distinct(Group),
    N_comparisons = n_distinct(`Comparison (group1/group2)`),
    .groups = "drop"
  ) %>%
  arrange(Genes, desc(N_sites))

print(poi_presence)
write.csv(poi_presence,
          "tables/step3_poi_presence_summary2.csv",
          row.names = FALSE)

# Summarise evidence for each PTM class 
# to create a formal selection justification table

selection_justification <- ptm_class_overview %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons) %>%
  group_by(PTM.ModificationTitle) %>%
  summarise(
    Mean_proteins_across_comparisons = round(
      mean(N_unique_proteins), 0),
    Max_mean_abs_log2FC = round(
      max(Mean_abs_log2FC), 3),
    N_key_comparisons_present = n(),
    Detected_on_POI = any(PTM.ModificationTitle %in% 
                            unique(ptm_poi$PTM.ModificationTitle)),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_mean_abs_log2FC))

print(selection_justification)
write.csv(selection_justification,
          "tables/step3_ptm_selection_justification.csv",
          row.names = FALSE)
