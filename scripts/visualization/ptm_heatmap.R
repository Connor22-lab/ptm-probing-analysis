#===============================================================================
# VISUAL ABSTRACT FIGURES - Simplified PTM Analysis Script
# Project: Computational discovery of novel PTMs in the inflammasome system
# Author: Connor Allison
# Date: January 2026
#===============================================================================

#===============================================================================
# SECTION 0: SETUP
#===============================================================================
rm(list = ls())
# install.packages(c("here", "tidyverse", "ggrepel", "patchwork"))

library(here)
library(tidyverse)
library(ggrepel)
library(patchwork)

# === PATHS ===
ptm_file <- here("data", "processed", "ptm_closed_search_corrected.tsv")
protein_file <- here("data", "processed", "ptm_closed_search_corrected_PG_candidates.tsv")
output_dir <- here("results", "figures", "closed_search")

dir.create(output_dir, showWarnings = FALSE)

#===============================================================================
# SECTION 1: LOAD DATA AND DEFINE POI LIST
#===============================================================================

#===============================================================================
# SECTION 1: LOAD DATA AND DEFINE POI LIST
#===============================================================================

ptm_data <- read.delim(ptm_file)
protein_data <- read.delim(protein_file)

cat("PTM data:", nrow(ptm_data), "rows\n")
cat("Protein data:", nrow(protein_data), "rows\n")

# Check Group column format
cat("\n=== Group column examples ===\n")
print(head(unique(ptm_data$Group), 10))

# Simple POI list
poi_list <- data.frame(
  Gene = c("NLRP3", "PYCARD", "CASP1", "GSDMD", "IL1B", "IL18",
           "IKBKG", "LYN", "TRIM25", "TRIM33", "CARD9","NEK7","BRCC3" ),
  UniProt = c("Q96P20", "Q9ULZ3", "P29466", "P57764", "P01584", "Q14116",
              "Q9Y6K9", "P07948", "Q14258", "Q9UPN9", "Q9H257", "Q8TDX7","P46736
")
)

#===============================================================================
# SECTION 2: FILTER AND NORMALIZE
#===============================================================================

artifact_ptms <- c("Carbamidomethyl (C)", "Oxidation (M)")

ptm_filtered <- ptm_data %>%
  filter(!PTM.ModificationTitle %in% artifact_ptms)

ptm_all <- ptm_data

cat("\nAfter removing artifacts:", nrow(ptm_filtered), "rows\n")

# Normalize to protein abundance
ptm_normalized <- ptm_filtered %>%
  left_join(
    protein_data %>%
      select(UniProtIds, Comparison..group1.group2., AVG.Log2.Ratio) %>%
      rename(Protein_Log2Ratio = AVG.Log2.Ratio),
    by = c("UniProtIds", "Comparison..group1.group2.")
  ) %>%
  mutate(Norm_Log2Ratio = AVG.Log2.Ratio - Protein_Log2Ratio)

cat("Normalization complete.\n")

#===============================================================================
# SECTION 3: PANEL 1A - DETECTION BAR CHART (by unique site count)
#===============================================================================

# Count UNIQUE sites per protein using Group column
# Group format: UniProt_ModificationResidue_...

# First, extract UniProt from Group column
site_counts_all <- ptm_all %>%
  filter(UniProtIds %in% poi_list$UniProt) %>%
  # Get unique sites (Group column)
  distinct(Group, UniProtIds) %>%
  # Count per protein
  group_by(UniProtIds) %>%
  summarise(Total_Sites = n(), .groups = "drop") %>%
  # Merge with gene names
  left_join(poi_list, by = c("UniProtIds" = "UniProt"))

# Ensure all POI are included (even with 0)
detection_df <- poi_list %>%
  left_join(site_counts_all %>% select(UniProtIds, Total_Sites), 
            by = c("UniProt" = "UniProtIds")) %>%
  mutate(Total_Sites = replace_na(Total_Sites, 0)) %>%
  arrange(desc(Total_Sites))

# Set factor order for plotting
detection_df$Gene <- factor(detection_df$Gene, levels = detection_df$Gene)

cat("\n=== Unique sites per protein (including artifacts) ===\n")
print(detection_df)

# Plot
p1a <- ggplot(detection_df, aes(x = Total_Sites, y = Gene)) +
  geom_col(fill = "#2C5F2D", width = 0.7) +
  geom_text(aes(label = Total_Sites), hjust = -0.3, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "PTM Sites Detected",
    subtitle = "Unique sites per protein (including artifacts)",
    x = "Number of Unique Sites",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Panel_1A_Detection.png"), p1a, 
       width = 6, height = 5, dpi = 300)
cat("\nPanel 1A saved.\n")

#===============================================================================
# SECTION 4: PANEL 1B - PTM COUNTS HEATMAP (Unique sites, relevant PTMs only)
#===============================================================================

# Count UNIQUE sites per protein per PTM type (excluding artifacts)
ptm_site_counts <- ptm_filtered %>%
  filter(UniProtIds %in% poi_list$UniProt) %>%
  # Get unique sites per PTM type
  distinct(Group, UniProtIds, PTM.ModificationTitle) %>%
  # Merge with gene names
  left_join(poi_list, by = c("UniProtIds" = "UniProt")) %>%
  # Count per protein per PTM type
  group_by(Gene, PTM.ModificationTitle) %>%
  summarise(Count = n(), .groups = "drop")

# Ensure all combinations exist
all_combos <- expand.grid(
  Gene = poi_list$Gene,
  PTM.ModificationTitle = unique(ptm_filtered$PTM.ModificationTitle)
)

ptm_counts_full <- all_combos %>%
  left_join(ptm_site_counts, by = c("Gene", "PTM.ModificationTitle")) %>%
  mutate(Count = replace_na(Count, 0))

# Order genes by total count
gene_order <- ptm_counts_full %>%
  group_by(Gene) %>%
  summarise(Total = sum(Count)) %>%
  arrange(desc(Total)) %>%
  pull(Gene)

ptm_counts_full$Gene <- factor(ptm_counts_full$Gene, levels = rev(gene_order))

cat("\n=== Unique PTM sites per protein (relevant PTMs only) ===\n")
print(ptm_counts_full %>% pivot_wider(names_from = PTM.ModificationTitle, values_from = Count))

# Plot
p1b <- ggplot(ptm_counts_full, aes(x = PTM.ModificationTitle, y = Gene, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#D94801", name = "Unique\nSites") +
  labs(
    title = "Unique PTM Sites",
    subtitle = "Biologically relevant PTMs only",
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "Panel_1B_PTM_Counts.png"), p1b, 
       width = 6, height = 5, dpi = 300)
cat("Panel 1B saved.\n")

#===============================================================================
# SECTION 5: PANEL 2 - VOLCANO PLOTS (P3C4 Model: Priming + Activation)
#===============================================================================

# Check available comparisons
cat("\n=== Available Comparisons ===\n")
comparisons <- unique(ptm_normalized$Comparison..group1.group2.)
print(comparisons)

# Find P3C4 comparisons
# Priming: P3C4 vs Unstim (or Unstim vs P3C4)
# Activation: P3C4_Nigericin vs P3C4 (or P3C4 vs P3C4_Nigericin)

# Define comparisons
priming_comp <- comparisons[grepl("P3C4.*Unstim|Unstim.*P3C4", comparisons, ignore.case = TRUE)][1]
activation_comp <- comparisons[grepl("P3C4.*Nigericin.*P3C4|P3C4.*P3C4.*Nigericin", comparisons, ignore.case = TRUE)][1]

cat("\nUsing:\n")
cat("Priming:", priming_comp, "\n")
cat("Activation:", activation_comp, "\n")

# Function to create improved volcano plot
make_volcano_v2 <- function(data, comparison, title) {
  
  plot_data <- data %>%
    filter(Comparison..group1.group2. == comparison) %>%
    filter(!is.na(Norm_Log2Ratio) & !is.na(Pvalue)) %>%
    mutate(
      neglog10P = -log10(Pvalue),
      log2FC = -Norm_Log2Ratio,  # Flip so positive = increase
      is_POI = UniProtIds %in% poi_list$UniProt
    ) %>%
    # Classify points
    mutate(
      Category = case_when(
        Pvalue < 0.05 & log2FC >= 0.58 ~ "Up",
        Pvalue < 0.05 & log2FC <= -0.58 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # POI data for labeling
  poi_data <- plot_data %>%
    filter(is_POI)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = log2FC, y = neglog10P)) +
    # All points
    geom_point(aes(color = Category), alpha = 0.5, size = 1.5) +
    
    # POI points (larger, purple, outlined)
    geom_point(data = poi_data, 
               color = "#9B59B6", size = 3, alpha = 0.9) +
    
    # Labels for POI only (purple)
    geom_text_repel(
      data = poi_data,
      aes(label = Genes),
      size = 2.8, color = "#9B59B6", fontface = "bold",
      max.overlaps = 20,
      segment.color = "#9B59B6",
      segment.size = 0.3,
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    
    # Threshold lines with labels
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    
    # Add text annotations for thresholds
    annotate("text", x = Inf, y = -log10(0.05), 
             label = "P = 0.05", 
             hjust = 1.1, vjust = -0.5, size = 2.5, color = "grey30") +
    annotate("text", x = 0.58, y = Inf, 
             label = "1.5-fold", 
             hjust = -0.1, vjust = 1.5, size = 2.5, color = "grey30", angle = 90) +
    annotate("text", x = -0.58, y = Inf, 
             label = "1.5-fold", 
             hjust = -0.1, vjust = -0.5, size = 2.5, color = "grey30", angle = 90) +
    
    # Colors: red = up, blue = down, grey = NS
    scale_color_manual(
      values = c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "grey70"),
      name = "Regulation",
      labels = c("Down" = "Down (P<0.05, FC<-1.5)", 
                 "NS" = "Not Significant", 
                 "Up" = "Up (P<0.05, FC>1.5)")
    ) +
    
    labs(
      title = title,
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10](P-value))
    ) +
    
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10)
    )
  
  return(p)
}

# Create both volcano plots
v1 <- make_volcano_v2(ptm_normalized, priming_comp, "Priming: Unstim → P3C4")
v2 <- make_volcano_v2(ptm_normalized, activation_comp, "Activation: P3C4 → P3C4+Nigericin")

# Combine with improved annotation
panel2 <- v1 + v2 +
  plot_annotation(
    title = "PTM Changes Across Inflammasome Activation",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
    )
  ) +
  plot_layout(guides = "collect") &  # Combine legends
  theme(legend.position = "bottom")

# Save
ggsave(file.path(output_dir, "Panel_2_Volcanos.png"), panel2, 
       width = 12, height = 6, dpi = 300)

cat("Panel 2 saved.\n")
cat("Changes made:\n")
cat("- Removed top 15 hit labels\n")
cat("- Only POI proteins are labeled (purple)\n")
cat("- Added threshold annotations on plot\n")
cat("- Improved legend with parameter details\n")
cat("- Updated subtitle with clearer information\n")

# SECTION 6: PANEL 3 - PTM DYNAMICS HEATMAP (CORRECTED)
#===============================================================================

# Get comparisons vs Unstim
unstim_comps <- comparisons[grepl("Unstim", comparisons, ignore.case = TRUE)]

# Filter POI data
dynamics_data <- ptm_normalized %>%
  filter(UniProtIds %in% poi_list$UniProt) %>%
  filter(Comparison..group1.group2. %in% unstim_comps) %>%
  filter(!is.na(Norm_Log2Ratio)) %>%
  left_join(poi_list, by = c("UniProtIds" = "UniProt")) %>%
  mutate(
    Condition = gsub(" / Unstim.*| / UNSTIM.*", "", Comparison..group1.group2.),
    Site = paste(Gene, PTM.ModificationTitle, sep = "\n")
  ) %>%
  # Average per unique site per condition
  group_by(Site, Gene, PTM.ModificationTitle, Condition) %>%
  summarise(
    # IMPORTANT: Flip the sign to match volcano plot interpretation
    Log2FC = -mean(Norm_Log2Ratio, na.rm = TRUE),  # ADDED NEGATIVE SIGN
    .groups = "drop"
  )

if(nrow(dynamics_data) > 0) {
  
  # Order conditions
  condition_order <- c("LPS", "P3C4", "Nigericin", "LPS_Nigericin", "P3C4_Nigericin")
  dynamics_data$Condition <- factor(dynamics_data$Condition, 
                                    levels = intersect(condition_order, unique(dynamics_data$Condition)))
  
  p3 <- ggplot(dynamics_data, aes(x = Condition, y = Site, fill = Log2FC)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#2166AC",      # Blue = downregulated in condition vs Unstim
      mid = "white", 
      high = "#B2182B",     # Red = upregulated in condition vs Unstim
      midpoint = 0, 
      name = "Log2FC\n(vs Unstim)",
      limits = c(-max(abs(dynamics_data$Log2FC), na.rm = TRUE), 
                 max(abs(dynamics_data$Log2FC), na.rm = TRUE))  # Symmetric scale
    ) +
    labs(
      title = "PTM Dynamics Across Conditions",
      x = "Condition", 
      y = "Protein - PTM"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, "Panel_3_Dynamics.png"), p3, 
         width = 10, height = 8, dpi = 300)
  cat("Panel 3 saved with CORRECTED Log2FC values (flipped to match volcano plots).\n")
  cat("Red = upregulated in condition vs Unstim\n")
  cat("Blue = downregulated in condition vs Unstim\n")
}

# VERIFICATION: Check TRIM33 in P3C4 condition
cat("\n=== VERIFICATION ===\n")
trim33_check <- dynamics_data %>%
  filter(grepl("TRIM33", Gene, ignore.case = TRUE)) %>%
  filter(Condition == "P3C4")

if(nrow(trim33_check) > 0) {
  cat("TRIM33 in P3C4 vs Unstim:\n")
  print(trim33_check %>% select(Gene, PTM.ModificationTitle, Condition, Log2FC))
  cat("\nPositive Log2FC means upregulated in P3C4 (matches volcano plot)\n")
} else 
  cat("TRIM33 not found in P3C4 condition\n")

#===============================================================================
# SECTION 7: SUMMARY
#===============================================================================

cat("\n========================================\n")
cat("DONE! Files saved to:", output_dir, "\n")
cat("- Panel_1A_Detection.png (bar chart by unique sites)\n")
cat("- Panel_1B_PTM_Counts.png (unique sites per PTM type)\n")
cat("- Panel_2_Volcanos.png (P3C4 priming + activation)\n")
cat("- Panel_3_Dynamics.png\n")
cat("========================================\n")







