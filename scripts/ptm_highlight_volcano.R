# ====================================================================
# PTM HIGHLIGHT VOLCANO (ONE PLOT + ONE CSV)
# ====================================================================
# Generates a volcano plot for ONE comparison, highlighting ALL sites
# of a chosen PTM modification. Also exports a CSV of all those sites.
# ====================================================================

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# -----------------------------
# USER SETTINGS (EDIT THESE)
# -----------------------------
setwd("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw")
ptm_file <- "PTM1.tsv"

# LOAD DATA
# -----------------------------
PTM_data <- as.data.frame(fread(ptm_file))
message("Data loaded: ", nrow(PTM_data), " PTM observations")

# -----------------------------
# Let's see what experimental comparisons are in the data
print("=== All comparisons in your dataset ===")
all_comparisons <- unique(PTM_data$`Comparison (group1/group2)`)
for(i in 1:length(all_comparisons)){
  print(paste(i, ":", all_comparisons[i]))
}
# -----------------------------
my_comparison <- "P3C4 / LPS"
my_ptm <- "GlyGly (K)"      # e.g., "Acetyl (K)", "Phospho (STY)" , "GlyGly (K)", "HexNAc (ST)", etc.
label_n <- 15               # number of highlighted points to label
out_dir <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/results" # where PNG + CSV will be saved

# thresholds (Spectronaut-style defaults)
fc_thr <- 0.585             # log2(1.5)
q_thr <- 0.05
# -----------------------------
# HELPER: SAFE FILENAME
# -----------------------------
safe_name <- function(x) gsub("/| |\\+|\\(|\\)|\\:|\\,", "_", x)

# -----------------------------
# MAIN FUNCTION
# -----------------------------
make_highlight_volcano <- function(data, comparison, ptm, out_dir=".", label_n=15,
                                  fc_thr=0.585, q_thr=0.05) {

  df <- data %>%
    filter(`Comparison (group1/group2)` == comparison) %>%
    mutate(
      Highlight = ifelse(PTM.ModificationTitle == ptm, ptm, "Other PTMs"),
      Significance = case_when(
        `AVG Log2 Ratio` >=  fc_thr & Qvalue <= q_thr ~ "Up",
        `AVG Log2 Ratio` <= -fc_thr & Qvalue <= q_thr ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  # Filter rows that can actually be plotted (prevents geom_point warning)
  df_plot <- df %>%
    filter(!is.na(Qvalue), Qvalue > 0, !is.na(`AVG Log2 Ratio`))

  # ---- CSV: ALL sites with this PTM (use df, not df_plot, so you keep full info) ----
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
      Qvalue,
      Pvalue,
      `# of Ratios`,
      Significance
    ) %>%
    arrange(Qvalue)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_csv <- file.path(
    out_dir,
    paste0("All_", safe_name(ptm), "_sites_", safe_name(comparison), ".csv")
  )

  write.csv(all_ptm_list, out_csv, row.names = FALSE)
  message("Saved CSV: ", out_csv, " (rows: ", nrow(all_ptm_list), ")")

  # ---- Select points to label (top by Qvalue within highlighted PTM) ----
  top_hits <- df_plot %>%
    filter(Highlight == ptm) %>%
    arrange(Qvalue) %>%
    slice_head(n = label_n)

  # ---- Volcano plot with highlighted PTM ----
  p <- ggplot(df_plot, aes(x = `AVG Log2 Ratio`, y = -log10(Qvalue))) +
    geom_point(aes(color = Highlight, size = Highlight, alpha = Highlight)) +
    scale_color_manual(values = c("Other PTMs" = "grey80", ptm = "darkgreen"), name = "PTM Type") +
    scale_size_manual(values = c("Other PTMs" = 0.5, ptm = 2)) +
    scale_alpha_manual(values = c("Other PTMs" = 0.3, ptm = 1)) +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(q_thr), linetype = "dashed") +
    labs(
      title = paste("All", ptm, "sites highlighted\n", comparison),
      x = "Log2 Fold Change",
      y = "-Log10 Q-value"
    ) +
    guides(size = "none", alpha = "none") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # Label directly from top_hits (NO join => no many-to-many warning)
  if (nrow(top_hits) > 0) {
    p <- p + geom_text_repel(
      data = top_hits,
      aes(label = Genes),
      size = 2.5,
      color = "darkgreen",
      max.overlaps = 30   # increase if you still see overlap warnings
    )
  }

  out_png <- file.path(
    out_dir,
    paste0("Volcano_Highlighted_", safe_name(ptm), "_", safe_name(comparison), ".png")
  )

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
  q_thr = q_thr
)

print(res$plot)
message("Done. Highlighted sites: ", res$n_sites)

