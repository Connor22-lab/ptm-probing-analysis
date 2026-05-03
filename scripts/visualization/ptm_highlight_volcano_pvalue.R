# ====================================================================
# PTM HIGHLIGHT VOLCANO (ONE PLOT + ONE CSV)
# ====================================================================
# Generates a volcano plot for ONE comparison, highlighting ALL sites
# of a chosen PTM modification. Also exports a CSV of all those sites.
# NOW: uses Pvalue for plotting & significance (not Qvalue)
# ====================================================================

rm(list=ls())

library(here)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# -----------------------------
# USER SETTINGS (EDIT THESE)
# -----------------------------
ptm_file <- here("data_raw", "ptm_closed_search.tsv")

my_comparison <- "Unstim / P3C4"
my_ptm <- "Acetyl (K)"      # e.g., "Acetyl (K)", "Phospho (STY)" , "GlyGly (K)", "HexNAc (ST)", etc.
label_n <- 15               # number of highlighted points to label

out_dir <- here("results", "figures", "closed_search")

# thresholds
fc_thr <- 0.585             # log2(1.5)
p_thr  <- 0.05              # use Pvalue threshold now (not q)

# -----------------------------
# LOAD DATA
# -----------------------------
PTM_data <- as.data.frame(fread(ptm_file))
message("Data loaded: ", nrow(PTM_data), " PTM observations")

# -----------------------------
# Let's see what experimental comparisons are in the data
# -----------------------------
print("=== All comparisons in your dataset ===")
all_comparisons <- unique(PTM_data$`Comparison (group1/group2)`)
for(i in 1:length(all_comparisons)){
  print(paste(i, ":", all_comparisons[i]))
}

# -----------------------------
# HELPER: SAFE FILENAME
# -----------------------------
safe_name <- function(x) gsub("[/\\s\\+\\(\\)\\:\\,]+", "_", x)

# -----------------------------
# MAIN FUNCTION
# -----------------------------
make_highlight_volcano <- function(data, comparison, ptm, out_dir=".", label_n=15,
                                   fc_thr=0.585, p_thr=0.05) {
  
  df <- data %>%
    filter(`Comparison (group1/group2)` == comparison) %>%
    mutate(
      Highlight = ifelse(PTM.ModificationTitle == ptm, ptm, "Other PTMs"),
      Significance = case_when(
        `AVG Log2 Ratio` >=  fc_thr & Pvalue <= p_thr ~ "Up",
        `AVG Log2 Ratio` <= -fc_thr & Pvalue <= p_thr ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  # Filter rows that can actually be plotted
  df_plot <- df %>%
    filter(!is.na(Pvalue), Pvalue > 0, !is.na(`AVG Log2 Ratio`))
  
  # ---- CSV: ALL sites with this PTM (keep full info incl Qvalue) ----
  all_ptm_list <- df %>%
    filter(Highlight == ptm) %>%
    select(
      `Comparison (group1/group2)`,
      Genes,
      UniProtIds,
      ProteinGroups,
      ProteinNames,
      ProteinDescriptions,
      Group,
      PTM.ModificationTitle,
      `AVG Log2 Ratio`,
      Pvalue,
      Qvalue,
      `# of Ratios`,
      Significance
    ) %>%
    arrange(Pvalue)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ---- Filenames: volcano_<PTMtype>_<Comparison>.png ----
  out_png <- file.path(
    out_dir,
    paste0("volcano_", safe_name(ptm), "_", safe_name(comparison), ".png")
  )
  
  out_csv <- file.path(
    out_dir,
    paste0("sites_", safe_name(ptm), "_", safe_name(comparison), ".csv")
  )
  
  write.csv(all_ptm_list, out_csv, row.names = FALSE)
  message("Saved CSV: ", out_csv, " (rows: ", nrow(all_ptm_list), ")")
  
  # ---- Select points to label (top by Pvalue within highlighted PTM) ----
  top_hits <- df_plot %>%
    filter(Highlight == ptm) %>%
    arrange(Pvalue) %>%
    slice_head(n = label_n)
  
  # ---- Volcano plot (Pvalue on Y) ----
  p <- ggplot(df_plot, aes(x = `AVG Log2 Ratio`, y = -log10(Pvalue))) +
    geom_point(aes(color = Highlight, size = Highlight, alpha = Highlight)) +
    scale_color_manual(values = c("Other PTMs" = "grey80", ptm = "darkgreen"), name = "PTM Type") +
    scale_size_manual(values = c("Other PTMs" = 0.5, ptm = 2)) +
    scale_alpha_manual(values = c("Other PTMs" = 0.3, ptm = 1)) +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_thr), linetype = "dashed") +
    labs(
      title = paste("All", ptm, "sites highlighted\n", comparison),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    guides(size = "none", alpha = "none") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (nrow(top_hits) > 0) {
    p <- p + geom_text_repel(
      data = top_hits,
      aes(label = Genes),
      size = 2.5,
      color = "darkgreen",
      max.overlaps = 30
    )
  }
  
  ggsave(out_png, plot = p, width = 22, height = 16, units = "cm", dpi = 300)
  message("Saved PNG: ", out_png)
  
  return(list(plot = p, csv = out_csv, png = out_png, n_sites = nrow(all_ptm_list)))
}

# -----------------------------
# RUN (ONE PLOT)
# -----------------------------
res <- make_highlight_volcano(
  data = PTM_data,
  comparison = my_comparison,
  ptm = my_ptm,
  out_dir = out_dir,
  label_n = label_n,
  fc_thr = fc_thr,
  p_thr = p_thr
)

print(res$plot)
message("Done. Highlighted sites: ", res$n_sites)
