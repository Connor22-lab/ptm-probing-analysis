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


# STEP 3: LOAD YOUR PTM DATA
# This is the same data file from Script 1
PTM_data <- as.data.frame(fread(here("data_raw", "PTM1.tsv")))

print(paste("Data loaded:", nrow(PTM_data), "PTM observations"))

PTM_data$`Comparison (group1/group2)` <-
  gsub("\\s+", " ", trimws(PTM_data$`Comparison (group1/group2)`))


# ---- PROJECT PATHS ----
PROJECT_ROOT <- dirname(getwd())
RUN_TAG <- "open_search"   # change to "open_search" when needed
OUT_DIR <- file.path(PROJECT_ROOT, "figures", RUN_TAG)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
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
  "Unstim / P3C4",              # Effect of LPS priming
  "P3C4 / P3C4+Nigericin")

P3C4_comparisons <- gsub("\\s+", " ", trimws(P3C4_comparisons))

print("\n=== We will analyze these 5 LPS-focused comparisons ===")
for(comp in P3C4_comparisons){
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
                              volcano_data$Pvalue <= 0.05] <- "Up"
  
  # Mark DOWN-regulated PTMs (decreased in numerator condition)
  volcano_data$Significance[volcano_data$`AVG Log2 Ratio` <= -0.585 & 
                              volcano_data$Pvalue <= 0.05] <- "Down"
  
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
    arrange(Pvalue) %>%                            # Sort by Q-value (lowest first)
    head(10)                                       # Take top 10
  
  # STEP 5: CREATE THE VOLCANO PLOT
  p <- ggplot(volcano_data, aes(x = `AVG Log2 Ratio`, y = -log10(Pvalue))) +
    
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
         y = "-Log10 P-value") +          # Y-axis label
    
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
  filename <- paste0("Volcano_", safe_name(comparison_name), "_", RUN_TAG, ".png")
  
  
  ggsave(filename = file.path(OUT_DIR, filename),               # Save with this filename
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
for(comp in P3C4_comparisons) {
  p <- make_volcano(comp)  # Create and save plot
  print(p)                 # Display plot in RStudio
}

print("\n=== PART A COMPLETE: 5 volcano plots saved ===")
print("Check your folder for Volcano_*.png files")

# ====================================================================
# PART C: HIGHLIGHT SPECIFIC PTM TYPES
# ====================================================================
# Now we want to look at specific modifications (like ubiquitination)
# and see ALL sites with that modification, highlighted on the volcano

print("\n========================================")
print("PART C: HIGHLIGHTING SPECIFIC PTM TYPES")
print("========================================")

# STEP 1: Choose which comparison to focus on
my_comparison <- "P3C4 / P3C4+Nigericin"# YOU CAN CHANGE THIS

print(paste("Focusing on:", my_comparison))

# STEP 2: See what PTM types are available
print("\n=== Available PTM types in your data ===")
available_ptms <- unique(PTM_data$PTM.ModificationTitle)

# Print them with numbers for easy reference
for(i in 1:length(available_ptms)){
  print(paste(i, ":", available_ptms[i]))
}

# STEP 3: Count PTMs in your chosen comparison
ptm_counts <- PTM_data %>%
  filter(`Comparison (group1/group2)` == my_comparison) %>%  # Filter to your comparison
  group_by(PTM.ModificationTitle) %>%                        # Group by PTM type
  summarise(
    Total = n(),                                             # Count total sites
    Significant = sum(Pvalue < 0.05 & abs(`AVG Log2 Ratio`) > 0.585)  # Count significant
  ) %>%
  arrange(desc(Total))  # Sort by most abundant

print("\n=== PTM counts in this comparison ===")
print(ptm_counts)

# STEP 4: Choose which PTM to highlight

my_ptm <- "Acetyl (K)"  # YOU CAN CHANGE THIS - use exact name from list above

print(paste("\nHighlighting:", my_ptm))

# STEP 5: Prepare data for highlighted plot
highlight_data <- PTM_data %>%
  filter(`Comparison (group1/group2)` == my_comparison)

# Create a new column that labels your PTM of interest vs others
highlight_data$Highlight <- ifelse(
  highlight_data$PTM.ModificationTitle == my_ptm,  # If PTM matches...
  my_ptm,                                          # ...label with PTM name
  "Other PTMs"                                     # ...otherwise "Other PTMs"
)

# Label significance separately
highlight_data$Significance <- "Not significant"
highlight_data$Significance[highlight_data$`AVG Log2 Ratio` >= 0.585 & 
                              highlight_data$Pvalue <= 0.05] <- "Up"
highlight_data$Significance[highlight_data$`AVG Log2 Ratio` <= -0.585 & 
                              highlight_data$Pvalue <= 0.05] <- "Down"

# STEP 6: Count highlighted PTMs
n_highlighted <- sum(highlight_data$Highlight == my_ptm)
n_sig_highlighted <- sum(highlight_data$Highlight == my_ptm & 
                           highlight_data$Significance != "Not significant")

print(paste("Total", my_ptm, "sites:", n_highlighted))
print(paste("Significant", my_ptm, "sites:", n_sig_highlighted))

# STEP 7: Get top sites for labeling (whether significant or not)
top_ptm_hits <- highlight_data %>%
  filter(Highlight == my_ptm) %>%  # Only your PTM
  arrange(Pvalue) %>%               # Sort by Q-value
  head(15)                          # Take top 15

# STEP 8: CREATE HIGHLIGHTED VOLCANO PLOT
p_highlight <- ggplot(highlight_data, 
                      aes(x = `AVG Log2 Ratio`, y = -log10(Pvalue))) +
  
  # Add points with different sizes and colors for highlighted PTM
  geom_point(aes(color = Highlight,    # Color by PTM type
                 size = Highlight,      # Size by PTM type
                 alpha = Highlight)) +  # Transparency by PTM type
  
  # Set colors - your PTM in green, others in grey
  scale_color_manual(values = c("Other PTMs" = "grey80", 
                                my_ptm = "darkgreen"),  # Green for your PTM
                     name = "PTM Type") +
  
  # Make your PTM bigger
  scale_size_manual(values = c("Other PTMs" = 0.5,   # Small dots for others
                               my_ptm = 2)) +          # Big dots for your PTM
  
  # Make your PTM more visible
  scale_alpha_manual(values = c("Other PTMs" = 0.3,  # Transparent others
                                my_ptm = 1)) +         # Solid your PTM
  
  # Add threshold lines
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # Labels
  labs(title = paste("All", my_ptm, "sites highlighted\n", my_comparison),
       x = "Log2 Fold Change",
       y = "-Log10 Q-value") +
  
  # Theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right") +
  
  # Hide size and alpha legends (keep only color legend)
  guides(size = "none", alpha = "none")

# STEP 9: Add labels for top sites
if(nrow(top_ptm_hits) > 0) {
  p_highlight <- p_highlight + 
    geom_text_repel(data = top_ptm_hits,        # Label top sites
                    aes(label = Genes),          # Use gene names
                    size = 2.5,                  # Text size
                    color = "darkgreen",         # Green text
                    max.overlaps = 20)           # Allow more labels
}

# Display the plot
print(p_highlight)

# STEP 10: Save highlighted plot
highlight_filename <- paste0(
  "Volcano_Highlighted_",
  safe_name(my_ptm), "_",
  safe_name(my_comparison), "_",
  RUN_TAG,
  ".png"
)


ggsave(filename = file.path(OUT_DIR, highlight_filename), 
       plot = p_highlight, 
       width = 22, 
       height = 16, 
       units = "cm", 
       dpi = 300)

print(paste("Saved:", highlight_filename))

# STEP 11: Save list of ALL sites with your PTM
all_ptm_list <- highlight_data %>%
  filter(Highlight == my_ptm) %>%  # Only your PTM
  select(Genes,                    # Select relevant columns
         ProteinGroups, 
         Group, 
         PTM.ModificationTitle, 
         `AVG Log2 Ratio`, 
         Qvalue, 
         Pvalue, 
         Significance) %>%
  arrange(Pvalue)                  # Sort by significance

# Save as CSV
csv_filename <- paste0(
  "All_",
  safe_name(my_ptm), "_sites_",
  safe_name(my_comparison), "_",
  RUN_TAG,
  ".csv"
)


write.csv(all_ptm_list, file.path(OUT_DIR, csv_filename), row.names = FALSE)

print(paste("Saved:", csv_filename))

# Show top 15 in console
print("\n=== Top 15 Sites (sorted by Q-value) ===")
print(head(all_ptm_list, 15))

