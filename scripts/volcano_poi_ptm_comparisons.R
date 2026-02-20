# ====================================================================
# SCRIPT 2: VOLCANO PLOTS FOR PTM ANALYSIS
# ====================================================================
# This script creates volcano plots to visualize PTM changes
# We'll focus on P3C4 (TLR1/2 agonist) and NLRP3 activation
# ====================================================================

# STEP 1: CLEAR WORKSPACE AND LOAD PACKAGES
rm(list=ls())  # Clear everything from previous sessions

library(data.table)  # For fast file reading
library(dplyr)       # For data manipulation
library(ggplot2)     # For making plots
library(ggrepel)     # For adding nice labels to plots


# STEP 2: SET WORKING DIRECTORY
# IMPORTANT: Change this to YOUR file location
setwd("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw")
out_dir <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/figures/open_search"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# STEP 3: LOAD YOUR PTM DATA
# This is the same data file from Script 1
PTM_data <- as.data.frame(fread("PTM1.tsv"))

print(paste("Data loaded:", nrow(PTM_data), "PTM observations"))

# -----------------------------
# LOAD POI LIST (TEXT FILE)
# -----------------------------
poi_file <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/docs/poi.txt"

POI <- readLines(poi_file, warn = FALSE) |> trimws()
POI <- POI[POI != ""]  # drop blank lines
message("Loaded POIs: ", length(POI))


# ====================================================================
# PART A: UNDERSTANDING YOUR COMPARISONS
# ====================================================================

# Let's see what experimental comparisons are in the data
print("=== All comparisons in your dataset ===")
all_comparisons <- unique(PTM_data$`Comparison (group1/group2)`)
for(i in 1:length(all_comparisons)){
  print(paste(i, ":", all_comparisons[i]))
}

# For this analysis, we focus on P3C4-related comparisons
# P3C4 = TLR1/2 agonist (priming signal)
# Nigericin = NLRP3 activator (activation signal)
# P3C4+Nigericin = Full inflammasome activation

P3C4_comparisons <- c(
"P3C4 / P3C4+Nigericin",      # P3C4 primed effect vs LPS primed under activation
"Unstim / P3C4"  # Unstim vs P3C4+LPS (priming effect
)

print("\n=== We will analyze these 5 LPS-focused comparisons ===")
for(comp in P3C4_comparisons){
  print(comp)
}

# ====================================================================
# PART B: CREATE VOLCANO PLOTS
# ====================================================================
has_any_poi <- function(x, poi_vec) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  sapply(strsplit(x, ";"), function(tokens) {
    any(trimws(tokens) %in% poi_vec)
  })
}

# Set the plot style - black and white theme looks professional
theme_set(theme_bw(base_size = 12))

# FUNCTION: make_volcano()
# This function creates one volcano plot for a given comparison
# INPUT: comparison_name (e.g., "Unstim / P3C4")
# OUTPUT: A volcano plot showing PTM changes

make_volcano <- function(comparison_name) {
  
  volcano_data <- PTM_data %>%
    filter(`Comparison (group1/group2)` == comparison_name)
  
  print(paste("\n--- Processing:", comparison_name, "---"))
  print(paste("Total PTM sites:", nrow(volcano_data)))
  
  # -----------------------------
  # SIGNIFICANCE CALLS (Q-value based, as in your original)
  # -----------------------------
  volcano_data$Significance <- "Not significant"
  
  volcano_data$Significance[volcano_data$`AVG Log2 Ratio` >= 0.585 &
                              volcano_data$Qvalue <= 0.05] <- "Up"
  
  volcano_data$Significance[volcano_data$`AVG Log2 Ratio` <= -0.585 &
                              volcano_data$Qvalue <= 0.05] <- "Down"
  
  n_up <- sum(volcano_data$Significance == "Up")
  n_down <- sum(volcano_data$Significance == "Down")
  n_not_sig <- sum(volcano_data$Significance == "Not significant")
  
  print(paste("Significant UP:", n_up))
  print(paste("Significant DOWN:", n_down))
  print(paste("Not significant:", n_not_sig))
  
  # -----------------------------
  # POI FLAGGING
  # Choose ONE (or combine) depending on what your POI file contains:
  # -----------------------------
  
  # Option A: POI file contains GENE symbols (matches the Genes column)
  # ---- POI flag ----
  volcano_data$IsPOI <- has_any_poi(volcano_data$Genes, POI)
  message("POI hits in this comparison: ", sum(volcano_data$IsPOI, na.rm = TRUE))
  
  # ---- Plot group: POI overrides significance for color ----
  volcano_data$PlotGroup <- volcano_data$Significance
  volcano_data$PlotGroup[volcano_data$IsPOI] <- "POI"
  
  # ---- label ONLY top 15 POIs by Qvalue ----
  poi_hits <- volcano_data %>%
    filter(IsPOI) %>%
    filter(!is.na(Qvalue), Qvalue > 0) %>%
    arrange(Qvalue) %>%
    slice_head(n = 15)
  
  # ---- Volcano plot ----
  p <- ggplot(volcano_data, aes(x = `AVG Log2 Ratio`, y = -log10(Qvalue))) +
    
    # First: plot NON-POIs (background)
    geom_point(
      data = subset(volcano_data, !IsPOI),
      aes(color = PlotGroup),
      size = 1,
      alpha = 0.65
    ) +
    
    # Second: plot POIs on top
    geom_point(
      data = subset(volcano_data, IsPOI),
      aes(color = PlotGroup),
      size = 1.3,
      alpha = 0.95
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "red",
        "Down" = "blue",
        "Not significant" = "grey",
        "POI" = "green3"
      ),
      breaks = c("Up", "Down", "POI", "Not significant")
    ) +
    
    geom_vline(xintercept = c(-0.585, 0.585),
               linetype = "dashed",
               color = "black") +
    
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "black") +
    
    labs(title = comparison_name,
         x = "Log2 Fold Change",
         y = "-Log10 Q-value") +
    
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "right")
  
  # ---- Add POI labels (top 15 only) ----
  if (nrow(poi_hits) > 0) {
    p <- p + geom_text_repel(
      data = poi_hits,
      aes(label = Genes),
      size = 2.7,
      max.overlaps = Inf,   # since you only have up to 15, this is safe
      box.padding = 0.4,
      segment.size = 0.2
    )
  }
  
filename <- file.path(out_dir, paste0("Volcano_POI_", gsub("/| |\\+", "_", comparison_name), ".png"))
  
  ggsave(
    filename,
    plot = p,
    width = 20,
    height = 15,
    units = "cm",
    dpi = 300
  )
  
  print(paste("Saved:", filename))
  
  return(p)
}

  
# ====================================================================
# STEP 8: GENERATE ALL 5 VOLCANO PLOTS
# ====================================================================

print("\n========================================")
print("CREATING VOLCANO PLOTS FOR P3C4 ANALYSIS")
print("========================================")

# Loop through each comparison and create a volcano plot
for(comp in P3C4_comparisons) {
  p <- make_volcano(comp)  # Create and save plot
  print(p)                 # Display plot in RStudio
}

print("\n=== PART A COMPLETE: 5 volcano plots saved ===")
print("Check your folder for Volcano_*.png files")
