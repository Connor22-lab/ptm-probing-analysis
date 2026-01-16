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

my_comparison <- "LPS / LPS+Nigericin"
my_ptm <- "HexNAc (ST)"      # e.g., "Acetyl (K)", "Phospho (STY)"
label_n <- 15               # number of highlighted points to label
out_dir <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/results" # where PNG + CSV will be saved

# thresholds (Spectronaut-style defaults)
fc_thr <- 0.585             # log2(1.5)
q_thr <- 0.05

# -----------------------------
# LOAD DATA
# -----------------------------
PTM_data <- as.data.frame(fread(ptm_file))
message("Data loaded: ", nrow(PTM_data), " PTM observations")

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

  # ---- CSV: ALL sites with this PTM ----
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

  out_csv <- file.path(
    out_dir,
    paste0("All_", safe_name(ptm), "_sites_", safe_name(comparison), ".csv")
  )

  write.csv(all_ptm_list, out_csv, row.names = FALSE)
  message("Saved CSV: ", out_csv, " (rows: ", nrow(all_ptm_list), ")")

  # ---- Select points to label (top by Qvalue within highlighted PTM) ----
  top_hits <- all_ptm_list %>% head(label_n)

  # ---- Volcano plot with highlighted PTM ----
  p <- ggplot(df, aes(x = `AVG Log2 Ratio`, y = -log10(Qvalue))) +
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

  if (nrow(top_hits) > 0) {
    # Need coordinates for labels; pull from df for those rows by matching key columns
    # (Genes + ProteinGroups + AVG Log2 Ratio + Qvalue is usually sufficient)
    label_df <- df %>%
      inner_join(
        top_hits %>% select(Genes, ProteinGroups, `AVG Log2 Ratio`, Qvalue),
        by = c("Genes", "ProteinGroups", "AVG Log2 Ratio", "Qvalue")
      )

    p <- p + geom_text_repel(
      data = label_df,
      aes(label = Genes),
      size = 2.5,
      color = "darkgreen",
      max.overlaps = 20
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

