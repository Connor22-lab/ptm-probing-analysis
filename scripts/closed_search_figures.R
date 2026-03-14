# ============================================================
# Closed Search Analysis Pipeline
# ============================================================
rm(list=ls())
library(tidyverse)

# Output folders
csv_out_closed <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/closed_search/csv"
fig_out_closed <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/closed_search/figures"

dir.create(csv_out_closed, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_out_closed, recursive = TRUE, showWarnings = FALSE)

# Load files
closed_raw <- read_tsv("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw/ptm_closed_search_noimputation.tsv")

protein_closed_raw <- read_tsv("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw/ptm_closed_search_noimputation_PG_candidates.tsv")

# Check
nrow(closed_raw)
colnames(closed_raw)
nrow(protein_closed_raw)

# Flip comparisons
closed_flipped <- closed_raw %>%
  mutate(
    `Comparison (group1/group2)` = paste(
      `Condition Denominator`, "/", `Condition Numerator`
    ),
    `AVG Log2 Ratio` = -`AVG Log2 Ratio`
  )
# Define key comparisons for closed search
key_comparisons_closed <- c(
  "LPS / Unstim",
  "LPS_Nigericin / LPS",
  "LPS_Nigericin / Unstim",
  "P3C4 / Unstim",
  "P3C4_Nigericin / P3C4",
  "P3C4_Nigericin / Unstim"
)
# Filter to key comparisons
closed_filtered <- closed_flipped %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons_closed)

nrow(closed_filtered)
# ============================================================
# DEFINE POI LIST
# ============================================================

poi_genes <- c(
  # Core inflammasome machinery
  "NLRP3", "PYCARD", "CASP1", "GSDMD", "IL1B", "IL18",
  # NF-kB priming
  "IKBKG", "IKBKB", "NFKB1", "MYD88", "IRAK1", "TRAF6",
  # Ubiquitin regulators
  "BRCC3", "TNFAIP3", "TRIM25", "TRIM33",
  # Positive findings
  "CARD9", "IFI16", "HSP90AB1", "VDAC1",
  # Direct NLRP3 activator
  "NEK7"
)

poi_list <- c(
  "Q96P20","Q9ULZ3","P29466","P57764","P01584","Q14116",
  "Q9Y6K9","O14920","P19838","Q99836","P51617","Q9Y4K3",
  "P46736","P21580","Q14258","Q9UPN9",
  "Q9BWT7","Q16666","P08238","P21796",
  "Q8TDX7"
)

# Refilter
closed_poi <- closed_filtered %>%
  filter(Genes %in% poi_genes | UniProtIds %in% poi_list)

# Confirm
closed_poi %>% distinct(Genes) %>% arrange(Genes)

# Check biological PTMs on all 49 POIs
closed_poi %>%
  filter(PTM.ModificationTitle %in% c(
    "Acetyl (K)", "GlyGly (K)", "Phospho (STY)"
  )) %>%
  count(Genes, PTM.ModificationTitle) %>%
  arrange(Genes)
# ============================================================
# FIGURE 3 PANEL A — All PTMs detected on POIs (incl. artefacts)
# ============================================================

library(tidyverse)

fig3a_data <- closed_poi %>%
  filter(Genes != "TNFAIP3") %>%
  distinct(Genes, PTM.ModificationTitle, Group) %>%
  mutate(
    Residue = str_extract(Group, "(?<=_)[A-Z][0-9]+"),
    Residue = case_when(
      Genes == "IKBKG" & Residue %in% c("K139", "K143") ~ "K139/K143",
      Genes == "TRIM25" & Residue %in% c("K392", "K402") ~ "K392/K402",
      TRUE ~ Residue
    ),
    Group_collapsed = paste0(Genes, "_", Residue)
  ) %>%
  distinct(Genes, PTM.ModificationTitle, Group_collapsed) %>%
  count(Genes, PTM.ModificationTitle) %>%
  mutate(PTM.ModificationTitle = factor(PTM.ModificationTitle,
                                        levels = c("Carbamidomethyl (C)", "Oxidation (M)",
                                                   "Phospho (STY)", "Acetyl (K)", "GlyGly (K)")))
fig3a <- ggplot(fig3a_data, aes(x = reorder(Genes, n),
                                y = n,
                                fill = PTM.ModificationTitle)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Carbamidomethyl (C)" = "#d3d3d3",
    "Oxidation (M)"       = "#a9a9a9",
    "Phospho (STY)"       = "#2166ac",
    "Acetyl (K)"          = "#d6604d",
    "GlyGly (K)"          = "#4dac26"
  )) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of unique modified sites",
    fill = "Modification",
    title = "Unique PTM sites detected on inflammasome POIs (all modifications)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

fig3a

ggsave("C:/Users/conno/Desktop/Final Year School/Final Year Project/closed_search/figures/fig3a_unique_sites_pois.png",
       fig3a, width = 10, height = 8, dpi = 300)

#Biological PTMs only
# Define biological category order
protein_order <- c(
  # Effectors
  "IL1B", "GSDMD",
  # NF-kB priming
  "IKBKG",
  # Ubiquitin regulators  
  "TRIM25", "TRIM33",
  # Kinases/regulators
  "CARD9", "IFI16", "HSP90AB1", "VDAC1"
)

# Updated data prep — collapses K139/K143
fig3b_data <- closed_poi %>%
  filter(PTM.ModificationTitle %in% c(
    "Acetyl (K)", "GlyGly (K)", "Phospho (STY)"
  )) %>%
  distinct(Genes, PTM.ModificationTitle, Group) %>%
  mutate(
    Residue = str_extract(Group, "(?<=_)[A-Z][0-9]+"),
    Residue = case_when(
      Genes == "IKBKG" & Residue %in% c("K139", "K143") ~ "K139/K143",
      Genes == "TRIM25" & Residue %in% c("K392", "K402") ~ "K392/K402",
      TRUE ~ Residue
    ), 
    Group_collapsed = paste0(Genes, "_", Residue)
  ) %>%
  distinct(Genes, PTM.ModificationTitle, Group_collapsed) %>%
  count(Genes, PTM.ModificationTitle)

# Complete grid — same as before
fig3b_complete <- expand.grid(
  Genes = protein_order,
  PTM.ModificationTitle = c("Acetyl (K)", "GlyGly (K)", "Phospho (STY)"),
  stringsAsFactors = FALSE
) %>%
  left_join(fig3b_data, by = c("Genes", "PTM.ModificationTitle")) %>%
  mutate(Genes = factor(Genes, levels = rev(protein_order)))
fig3b <- ggplot(fig3b_complete, aes(x = PTM.ModificationTitle,
                                    y = Genes,
                                    fill = n)) +
  geom_tile(colour = "black", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(n), "", n)),
            colour = ifelse(fig3b_complete$n == 1, "black", "white"),
            size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#f4a582", high = "#b2182b",
                      na.value = "white",
                      breaks = c(1, 2),
                      labels = c("1", "2"),
                      limits = c(1, 2)) +
  labs(
    x = "Modification type",
    y = NULL,
    fill = "Unique sites",
    title = "Biologically relevant PTM sites on inflammasome POIs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )

fig3b

ggsave("C:/Users/conno/Desktop/Final Year School/Final Year Project/closed_search/figures/fig3b_biological_ptms_tile.png",
       fig3b, width = 7, height = 5, dpi = 300)


# Load protein candidates 

protein_closed <- read_tsv("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw/ptm_closed_search_noimputation_PG_candidates.tsv")

# Check structure

colnames(protein_closed)

nrow(protein_closed)

# Check comparison format matches PTM file

protein_closed %>% 
  
  distinct(`Comparison (group1/group2)`) %>% 
  
  print(n = 20)
# Flip comparisons
protein_closed <- protein_closed %>%
  mutate(
    `Comparison (group1/group2)` = str_replace(
      `Comparison (group1/group2)`,
      "^(.*) / (.*)$", "\\2 / \\1"
    ),
    `AVG Log2 Ratio` = `AVG Log2 Ratio` * -1
  )

# Filter to key comparisons
key_comparisons_closed <- c(
  "LPS / Unstim",
  "LPS_Nigericin / LPS",
  "LPS_Nigericin / Unstim",
  "P3C4 / Unstim",
  "P3C4_Nigericin / P3C4",
  "P3C4_Nigericin / Unstim"
)

protein_closed_filtered <- protein_closed %>%
  filter(`Comparison (group1/group2)` %in% key_comparisons_closed)

# Check
nrow(protein_closed_filtered)
protein_closed_filtered %>% distinct(`Comparison (group1/group2)`)
# ============================================================
# PROTEIN ABUNDANCE CORRECTION — CLOSED SEARCH
# ============================================================

# Extract protein log2FC lookup
protein_lookup <- protein_closed_filtered %>%
  select(UniProtIds, `Comparison (group1/group2)`, 
         protein_log2fc = `AVG Log2 Ratio`) %>%
  distinct()

# Join protein FC to PTM data and correct
closed_corrected <- closed_filtered %>%
  left_join(protein_lookup, 
            by = c("UniProtIds", "Comparison (group1/group2)")) %>%
  mutate(corrected_log2fc = `AVG Log2 Ratio` - protein_log2fc)

# Check match rate
n_total <- nrow(closed_corrected)
n_matched <- closed_corrected %>% filter(!is.na(corrected_log2fc)) %>% nrow()
n_unmatched <- closed_corrected %>% filter(is.na(corrected_log2fc)) %>% nrow()

cat("Total rows:", n_total, "\n")
cat("Matched:", n_matched, "\n")
cat("Unmatched:", n_unmatched, "\n")
cat("Match rate:", round(n_matched/n_total * 100, 1), "%\n")

# ============================================================
# REFILTER TO POIS WITH CORRECTED LOG2FC
# ============================================================

closed_poi <- closed_corrected %>%
  filter(Genes %in% poi_genes | UniProtIds %in% poi_list)

# Extract biological PTMs for dynamics heatmap
fig4a_data <- closed_poi %>%
  filter(PTM.ModificationTitle %in% c(
    "Acetyl (K)", "GlyGly (K)", "Phospho (STY)"
  )) %>%
  distinct(Genes, Group, PTM.ModificationTitle,
           `Comparison (group1/group2)`, corrected_log2fc, Pvalue) %>%
  mutate(
    # Clean site label
    Site = paste0(Genes, "_", str_extract(Group, "(?<=_)[A-Z][0-9]+")),
    # Significance flag
    sig = ifelse(Pvalue < 0.05, "*", "")
  )

# Check structure
fig4a_data %>%
  select(Site, PTM.ModificationTitle, 
         `Comparison (group1/group2)`, 
         corrected_log2fc, sig) %>%
  print(n = 30)
# ============================================================
# FIGURE 4 PANEL A — UPDATED DYNAMICS HEATMAP
# Grey missing cells + K139/K143 collapsed
# ============================================================

# Define comparison order
comp_order <- c(
  "LPS / Unstim",
  "LPS_Nigericin / LPS",
  "LPS_Nigericin / Unstim",
  "P3C4 / Unstim",
  "P3C4_Nigericin / P3C4",
  "P3C4_Nigericin / Unstim"
)

# Extract biological PTMs and collapse sites
fig4a_data <- closed_poi %>%
  filter(PTM.ModificationTitle %in% c(
    "Acetyl (K)", "GlyGly (K)", "Phospho (STY)"
  )) %>%
  distinct(Genes, Group, PTM.ModificationTitle,
           `Comparison (group1/group2)`, corrected_log2fc, Pvalue) %>%
  mutate(
    Residue = str_extract(Group, "(?<=_)[A-Z][0-9]+"),
    Residue = case_when(
      Genes == "IKBKG"  & Residue %in% c("K139", "K143") ~ "K139/K143",
      Genes == "TRIM25" & Residue %in% c("K392", "K402") ~ "K392/K402",
      TRUE ~ Residue
    ),
    mod_short = case_when(
      PTM.ModificationTitle == "Acetyl (K)"    ~ "Ac",
      PTM.ModificationTitle == "GlyGly (K)"   ~ "Ub",
      PTM.ModificationTitle == "Phospho (STY)" ~ "Ph"
    ),
    Site = paste0(Genes, "_", Residue, " (", mod_short, ")"),
    sig  = ifelse(Pvalue < 0.05, "*", "")
  ) %>%
  distinct(Site, PTM.ModificationTitle,
           `Comparison (group1/group2)`, corrected_log2fc, sig)
# Create complete grid
all_sites <- unique(fig4a_data$Site)

complete_grid <- expand.grid(
  Site = all_sites,
  `Comparison (group1/group2)` = comp_order,
  stringsAsFactors = FALSE
)
site_order <- c(
  "IKBKG_K139/K143 (Ub)",
  "TRIM25_K392/K402 (Ac)", "TRIM33_K953 (Ac)",
  "CARD9_K60 (Ac)",
  "IL1B_K171 (Ac)",
  "IFI16_S153 (Ph)", "HSP90AB1_S255 (Ph)",
  "VDAC1_K224 (Ac)", "VDAC1_K256 (Ub)"
)
# Join — missing cells become NA = grey
fig4a_plot <- complete_grid %>%
  left_join(fig4a_data, by = c("Site", "Comparison (group1/group2)")) %>%
  mutate(
    `Comparison (group1/group2)` = factor(
      `Comparison (group1/group2)`, levels = comp_order),
    Site = factor(Site, levels = rev(site_order))
  )
# Plot
fig4a <- ggplot(fig4a_plot, aes(
  x = `Comparison (group1/group2)`,
  y = Site,
  fill = corrected_log2fc)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(sig), "", sig)), 
            colour = "black", size = 4) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    na.value = "grey85",        # grey for not detected
    name = "Corrected\nLog2FC"
  ) +
  labs(
    x = NULL,
    y = "PTM site",
    title = "PTM dynamics across inflammasome activation conditions"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.caption = element_text(size = 8, hjust = 0)
  )

fig4a

ggsave("C:/Users/conno/Desktop/Final Year School/Final Year Project/closed_search/figures/fig4a_dynamics_heatmap.png",
       fig4a, width = 9, height = 6, dpi = 300)
# ============================================================
# FIGURE 4 PANEL B — PTM Site Annotation Table
# ============================================================

library(tidyverse)
library(gridExtra)
library(grid)

# Table data
fig4b_data <- data.frame(
  Protein = c("IL1B", "IKBKG", "VDAC1", "CARD9"),
  PTM_Type = c("Acetylation", "Ubiquitination", "Ubiquitination", "Acetylation"),
  Site = c("K171", "K139/K143", "K256", "K60"),
  Localisation_Score = c("1", "< 0.5", "1", "1"),
  Database_Status = c("Not annotated", "Annotated", "Not annotated", "Not annotated"),
  stringsAsFactors = FALSE
)

# Colour scheme matching your existing table
header_fill  <- "#E07B2A"   # orange header
row_fill_odd <- "#F9DEC9"   # light salmon rows
row_fill_even<- "#F5CBA7"   # slightly darker alternating

# Build table as grob
table_grob <- tableGrob(
  fig4b_data,
  rows = NULL,
  cols = c("Protein", "PTM Type", "Site", "Localisation\nScore", "Database\nStatus"),
  theme = ttheme_minimal(
    core = list(
      bg_params = list(
        fill = c(row_fill_odd, row_fill_even),
        col  = "white",
        lwd  = 1.5
      ),
      fg_params = list(
        fontsize  = 12,
        fontface  = "plain",
        col       = "black"
      )
    ),
    colhead = list(
      bg_params = list(
        fill = header_fill,
        col  = "white",
        lwd  = 1.5
      ),
      fg_params = list(
        fontsize  = 12,
        fontface  = "bold",
        col       = "white"
      )
    )
  )
)

# Save
png("C:/Users/conno/Desktop/Final Year School/Final Year Project/closed_search/figures/fig4b_site_annotation_table.png",
    width = 2800, height = 700, res = 300)
grid.draw(table_grob)
dev.off()
