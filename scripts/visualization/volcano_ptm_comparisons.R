# ====================================================================
# SCRIPT 2: VOLCANO PLOTS FOR PTM ANALYSIS
# ====================================================================
# This script creates volcano plots to visualize PTM changes
# We'll focus on P3C4 (TLR1/2 agonist) and NLRP3 activation
# ====================================================================

# STEP 1: CLEAR WORKSPACE AND LOAD PACKAGES
rm(list=ls())  # Clear everything from previous sessions

library(here)
library(data.table)  # For fast file reading
library(dplyr)       # For data manipulation
library(ggplot2)     # For making plots
library(ggrepel)     # For adding nice labels to plots


# STEP 2: SET WORKING DIRECTORY
# IMPORTANT: Change this to YOUR file location


# STEP 3: LOAD YOUR PTM DATA
# This is the same data file from Script 1
PTM_data <- as.data.frame(fread(here("data_raw", "PTM1.tsv")))

print(paste("Data loaded:", nrow(PTM_data), "PTM observations"))

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

LPS+P3C4_comparisons <- c(
  "P3C4+Nigericin / LPS+Nigericin",              # Difference in activation depending on priming
  "P3C4+Nigericin / LPS",           # LPS primed vs P3C4 primed under activation
  "P3C4 / LPS",    # Difference in priming without activation
  "P3C4 / LPS+Nigericin",      # P3C4 primed effect vs LPS primed under activation
  "Unstim / P3C4+LPS"  # Unstim vs P3C4+LPS (priming effect
)

print("\n=== We will analyze these 5 LPS-focused comparisons ===")
for(comp in LPS+P3C4_comparisons){
  print(comp)
}

# ====================================================================
# PART B: CREATE VOLCANO PLOTS
# ====================================================================

# Set the plot style - black and white theme looks professional
theme_set(theme_bw(base_size = 12))

# FUNCTION: make_volcano()
# This function creates one volcano plot for a given comparison
# INPUT: comparison_name (e.g., "Unstim / P3C4")
# OUTPUT: A volcano plot showing PTM changes

make_volcano <- function(comparison_name) {
  
  # STEP 1: Filter data for this specific comparison
  volcano_data <- PTM_data %>%
    filter(`Comparison (group1/group2)` == comparison_name)
  
  print(paste("\n--- Processing:", comparison_name, "---"))
  print(paste("Total PTM sites:", nrow(volcano_data)))
  
  # STEP 2: Label each PTM as significant or not
  # Spectronaut defaults: Log2FC > 0.585 (= 1.5x fold change), Q-value < 0.05
  
  volcano_data$Significance <- "Not significant"  # Start with all "not significant"
  
  # Mark UP-regulated PTMs (increased in numerator condition)
  volcano_data$Significance[volcano_data$`AVG Log2 Ratio` >= 0.585 & 
                              volcano_data$Qvalue <= 0.05] <- "Up"
  
  # Mark DOWN-regulated PTMs (decreased in numerator condition)
  volcano_data$Significance[volcano_data$`AVG Log2 Ratio` <= -0.585 & 
                              volcano_data$Qvalue <= 0.05] <- "Down"
  
  # STEP 3: Count how many PTMs are significantly changed
  n_up <- sum(volcano_data$Significance == "Up")
  n_down <- sum(volcano_data$Significance == "Down")
  n_not_sig <- sum(volcano_data$Significance == "Not significant")
  
  print(paste("Significant UP:", n_up))
  print(paste("Significant DOWN:", n_down))
  print(paste("Not significant:", n_not_sig))
  
  # STEP 4: Get the top 10 most significant PTMs for labeling
  top_hits <- volcano_data %>%
    filter(Significance != "Not significant") %>%  # Only significant ones
    arrange(Qvalue) %>%                            # Sort by Q-value (lowest first)
    head(10)                                       # Take top 10
  
  # STEP 5: CREATE THE VOLCANO PLOT
  p <- ggplot(volcano_data, aes(x = `AVG Log2 Ratio`, y = -log10(Qvalue))) +
    
    # Add points for each PTM site
    geom_point(aes(color = Significance),  # Color by significance
               size = 1,                   # Size of dots
               alpha = 0.6) +              # Transparency (0.6 = 60% opaque)
    
    # Set colors: red = up, blue = down, grey = not significant
    scale_color_manual(values = c("Up" = "red", 
                                  "Down" = "blue", 
                                  "Not significant" = "grey")) +
    
    # Add vertical lines at fold change thresholds (±0.585)
    geom_vline(xintercept = c(-0.585, 0.585), 
               linetype = "dashed",        # Dashed line
               color = "black") +
    
    # Add horizontal line at significance threshold (Q = 0.05)
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed", 
               color = "black") +
    
    # Add labels
    labs(title = comparison_name,          # Title at top
         x = "Log2 Fold Change",          # X-axis label
         y = "-Log10 Q-value") +          # Y-axis label
    
    # Customize appearance
    theme(plot.title = element_text(hjust = 0.5,     # Center title
                                    face = "bold",    # Bold title
                                    size = 11),       # Font size
          legend.position = "right")                 # Legend on right
  
  # STEP 6: Add gene names for top significant PTMs
  if(nrow(top_hits) > 0) {  # Only if there are significant hits
    p <- p + geom_text_repel(
      data = top_hits,                # Use top_hits data
      aes(label = Genes),             # Label with gene names
      size = 2.5,                     # Text size
      max.overlaps = 10,              # Prevent too many overlapping labels
      box.padding = 0.5,              # Space around text
      segment.size = 0.2)             # Line connecting label to point
  }
  
  # STEP 7: Save the plot
  # Create filename by replacing special characters with underscores
  filename <- paste0("Volcano_", gsub("/| |\\+", "_", comparison_name), ".png")
  
  ggsave(filename,              # Save with this filename
         plot = p,              # Save this plot
         width = 20,            # 20 cm wide
         height = 15,           # 15 cm tall
         units = "cm",          # Units in centimeters
         dpi = 300)             # High resolution (300 dots per inch)
  
  print(paste("Saved:", filename))
  
  # Return the plot so we can view it
  return(p)
}

# ====================================================================
# STEP 8: GENERATE ALL 5 VOLCANO PLOTS
# ====================================================================

print("\n========================================")
print("CREATING VOLCANO PLOTS FOR P3C4 ANALYSIS")
print("========================================")

# Loop through each comparison and create a volcano plot
for(comp in LPS+P3C4_comparisons) {
  p <- make_volcano(comp)  # Create and save plot
  print(p)                 # Display plot in RStudio
}

print("\n=== PART A COMPLETE: 5 volcano plots saved ===")
print("Check your folder for Volcano_*.png files")