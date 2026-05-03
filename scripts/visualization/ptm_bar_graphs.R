# ====================================================================
# SCRIPT 1: LOADING DATA AND EXPLORING PTM TYPES
# ====================================================================
# This script loads your PTM data and creates three basic bar graphs
# to understand what post-translational modifications (PTMs) were detected
# ====================================================================

# STEP 1: CLEAR YOUR WORKSPACE
# This removes any old data/variables from previous sessions
rm(list=ls())

# STEP 2: LOAD REQUIRED PACKAGES
# These are like toolboxes that give R extra functions
library(data.table)  # For fast reading of large files
library(ggplot2)     # For making nice graphs
library(dplyr)       # For data manipulation (filtering, summarizing, etc.)

# STEP 3: SET YOUR WORKING DIRECTORY
# This tells R where to find your files and save your outputs
# IMPORTANT: Change this path to match YOUR computer


# STEP 4: LOAD YOUR PTM DATA
# fread() is a fast function to read tab-separated files (.tsv)
# We convert it to a data frame (a table format R understands)
ptm_data <- as.data.frame(fread(here("data_raw", "PTM1.tsv")))

# STEP 5: LOOK AT YOUR DATA
# This shows you the first 6 rows of your data
print("First few rows of data:")
head(ptm_data)

# This tells you how many rows (observations) and columns (variables) you have
print("Dataset dimensions:")
dim(ptm_data)  # Should show 828005 rows, 21 columns

# This lists all column names
print("Column names:")
colnames(ptm_data)

# STEP 6: SUMMARIZE PTM TYPES
# We want to count how many times each PTM type appears
# group_by() groups data by PTM type
# summarise() calculates statistics for each group
# n() counts the number of rows in each group
# n_distinct() counts unique values

ptm_summary <- ptm_data %>%                           # Take ptm_data, then...
  group_by(PTM.ModificationTitle) %>%                 # Group by PTM type, then...
  summarise(                                          # Calculate these for each group:
    Count = n(),                                      # Total number of observations
    Unique_Proteins = n_distinct(ProteinGroups),      # Number of different proteins
    Unique_Sites = n_distinct(Group)                  # Number of different modification sites
  ) %>%
  arrange(desc(Count))                                # Sort by Count (highest first)

# Look at the summary
print("PTM Summary:")
print(ptm_summary)

# Save this summary as a CSV file for your records
write.csv(ptm_summary, "PTM_Summary_Table.csv", row.names = FALSE)

# ====================================================================
# STEP 7: CREATE BAR GRAPH 1 - TOTAL OBSERVATIONS PER PTM
# ====================================================================

# Set up the plot theme (makes graphs look nice and professional)
theme_set(theme_bw(base_size = 12))  # Black and white theme, font size 12

# Create the plot
graph1 <- ggplot(ptm_summary,                                 # Use ptm_summary data
                 aes(x = reorder(PTM.ModificationTitle, -Count),  # x-axis: PTMs ordered by count
                     y = Count)) +                            # y-axis: Count values
  geom_bar(stat = "identity",                                 # Make bar chart (height = Count)
           fill = "steelblue",                                # Color bars blue
           alpha = 0.8) +                                     # Make slightly transparent
  geom_text(aes(label = Count),                               # Add numbers on top of bars
            vjust = -0.5,                                     # Position slightly above bars
            size = 3) +                                       # Text size
  labs(title = "Total PTM Observations by Type",              # Main title
       x = "PTM Type",                                        # x-axis label
       y = "Number of Observations") +                        # y-axis label
  theme(axis.text.x = element_text(angle = 90,                # Rotate x-axis labels 90 degrees
                                   vjust = 0.5,               # Center them
                                   hjust = 1,                 # Align to right
                                   size = 9),                 # Font size
        plot.title = element_text(hjust = 0.5,                # Center the title
                                  face = "bold")) +            # Make title bold
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))   # Add space above bars for numbers

# Display the plot
print(graph1)

# Save the plot as a high-quality image
ggsave("BarGraph_1_Total_Observations.png",    # Filename
       plot = graph1,                           # Which plot to save
       width = 30,                              # Width in cm
       height = 18,                             # Height in cm
       units = "cm",                            # Units
       dpi = 300)                               # Resolution (300 = publication quality)

# ====================================================================
# STEP 8: CREATE BAR GRAPH 2 - UNIQUE SITES PER PTM
# ====================================================================

graph2 <- ggplot(ptm_summary,                                       # Use ptm_summary data
                 aes(x = reorder(PTM.ModificationTitle, -Unique_Sites),  # x: PTMs ordered by sites
                     y = Unique_Sites)) +                           # y: Unique sites
  geom_bar(stat = "identity",                                       # Bar chart
           fill = "coral",                                          # Orange color
           alpha = 0.8) +                                           # Transparency
  geom_text(aes(label = Unique_Sites),                              # Add numbers
            vjust = -0.5,                                           # Above bars
            size = 3) +                                             # Text size
  labs(title = "Unique Modified Sites by PTM Type",                 # Title
       x = "PTM Type",                                              # x label
       y = "Number of Unique Sites") +                              # y label
  theme(axis.text.x = element_text(angle = 90,                      # Rotate labels
                                   vjust = 0.5, 
                                   hjust = 1,
                                   size = 9),
        plot.title = element_text(hjust = 0.5,                      # Center title
                                  face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))         # Space for numbers

# Display
print(graph2)

# Save
ggsave("BarGraph_2_Unique_Sites.png", 
       plot = graph2, 
       width = 30, 
       height = 18, 
       units = "cm", 
       dpi = 300)

# ====================================================================
# STEP 9: CREATE BAR GRAPH 3 - UNIQUE PROTEINS PER PTM
# ====================================================================

graph3 <- ggplot(ptm_summary,                                        # Use ptm_summary data
                 aes(x = reorder(PTM.ModificationTitle, -Unique_Proteins),  # x: PTMs by proteins
                     y = Unique_Proteins)) +                         # y: Unique proteins
  geom_bar(stat = "identity",                                        # Bar chart
           fill = "darkseagreen",                                    # Green color
           alpha = 0.8) +                                            # Transparency
  geom_text(aes(label = Unique_Proteins),                            # Add numbers
            vjust = -0.5,                                            # Above bars
            size = 3) +                                              # Text size
  labs(title = "Unique Proteins Modified by PTM Type",               # Title
       x = "PTM Type",                                               # x label
       y = "Number of Unique Proteins") +                            # y label
  theme(axis.text.x = element_text(angle = 90,                       # Rotate labels
                                   vjust = 0.5, 
                                   hjust = 1,
                                   size = 9),
        plot.title = element_text(hjust = 0.5,                       # Center title
                                  face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))          # Space for numbers

# Display
print(graph3)

# Save
ggsave("BarGraph_3_Unique_Proteins.png", 
       plot = graph3, 
       width = 30, 
       height = 18, 
       units = "cm", 
       dpi = 300)

# ====================================================================
# STEP 10: PRINT SUMMARY TO CONSOLE
# ====================================================================

print("========================================")
print("SUMMARY OF YOUR PTM DATA")
print("========================================")
print(paste("Total PTM types detected:", nrow(ptm_summary)))
print(paste("Total observations:", sum(ptm_summary$Count)))
print("")
print("Top 5 most abundant PTMs:")
print(head(ptm_summary, 5))
print("========================================")
print("SCRIPT COMPLETE!")
print("Check your folder for 3 bar graphs and 1 CSV file")
print("========================================")

# ====================================================================
# END OF SCRIPT 1
# ====================================================================
# You should now have:
# - PTM_Summary_Table.csv (table with all PTM counts)
# - BarGraph_1_Total_Observations.png
# - BarGraph_2_Unique_Sites.png  
# - BarGraph_3_Unique_Proteins.png
# ====================================================================
