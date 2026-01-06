# =========================
# PTM1 Volcano Pipeline
# =========================

# 0) Working directory (recommended: point to project root, not data_raw)
# If you insist on data_raw, keep as is.
setwd("C:/Users/conno/Desktop/Final Year School/Final Year Project")

# 1) Packages (install missing, then load)
pkgs <- c("data.table","dplyr","ggplot2","ggrepel","stringr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# 2) Parameters
log2fc_cut <- 0.26
q_cut <- 0.10

# Priority PTMs (confirmed in your dataset)
ptm_keep <- c("Phospho (STY)", "GlyGly (K)", "Acetyl (K)", "HexNAc (ST)")

# Comparisons to start with (will auto-drop any not present)
comparisons_keep <- c(
  "Unstim / LPS",
  "LPS / LPS+Nigericin",
  "Unstim / LPS+Nigericin",
  "Unstim / P3C4",
  "P3C4 / P3C4+Nigericin",
  "Unstim / P3C4+Nigericin"
)
# 3) I/O locations (keeps outputs organised)
data_file <- file.path("data_raw", "PTM1.tsv")
poi_file  <- file.path("data_raw", "POI.csv")  # move POI.csv here OR change path
out_fig   <- file.path("results", "figures")
out_tab   <- file.path("results", "tables")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

# 4) Load data + quick triage
df <- data.table::fread(data_file)

message("Rows x Cols: ", nrow(df), " x ", ncol(df))
message("PTMs present: ", length(unique(df$PTM.ModificationTitle)))
message("Comparisons present: ", length(unique(df$`Comparison (group1/group2)`)))

# Ensure requested PTMs and comparisons exist (prevents typos)
ptm_keep <- intersect(ptm_keep, unique(df$PTM.ModificationTitle))
comparisons_keep <- intersect(comparisons_keep, unique(df$`Comparison (group1/group2)`))

if (length(ptm_keep) == 0) stop("None of the requested PTMs were found in PTM.ModificationTitle.")
if (length(comparisons_keep) == 0) stop("None of the requested comparisons were found in Comparison (group1/group2).")

# 5) Load POI list (UniProt accessions)
poi <- data.table::fread(poi_file)
if (!"Entry" %in% names(poi)) stop("POI.csv must contain a column named 'Entry' with UniProt accessions.")
poi_list <- unique(poi$Entry)

# Robust matching for semicolon-delimited UniProtIds
poi_regex <- paste0("(^|;)", paste(poi_list, collapse="|"), "(;|$)")

# 6) Prep function: creates volcano variables + POI tag
prep_for_volcano <- function(d) {
  d %>%
    mutate(
      Qvalue = as.numeric(Qvalue),
      `AVG Log2 Ratio` = as.numeric(`AVG Log2 Ratio`),
      Qvalue_safe = pmax(Qvalue, 1e-300),
      neglog10Q = -log10(Qvalue_safe),
      Enrichment = case_when(
        `AVG Log2 Ratio` >=  log2fc_cut & Qvalue < q_cut ~ "UP",
        `AVG Log2 Ratio` <= -log2fc_cut & Qvalue < q_cut ~ "DOWN",
        TRUE ~ "NO"
      ),
      Genes = if_else(is.na(Genes) | Genes == "", as.character(ProteinGroups), as.character(Genes)),
      Label = if_else(stringr::str_detect(UniProtIds, poi_regex), "poi", "Other")
    ) %>%
    filter(!is.na(`AVG Log2 Ratio`), !is.na(Qvalue))
}

# 7) Plot function
make_volcano <- function(d, title, label_poi_only_sig = TRUE) {
  
  label_data <- if (label_poi_only_sig) {
    d %>% filter(Label == "poi", Enrichment != "NO")
  } else {
    d %>% filter(Label == "poi")
  }
  
  ggplot(d, aes(x = `AVG Log2 Ratio`, y = neglog10Q)) +
    geom_vline(
      xintercept = c(-log2fc_cut, 0, log2fc_cut),
      linetype = "dashed",
      color = c("black","darkgrey","black")
    ) +
    geom_hline(yintercept = -log10(q_cut), linetype = "dashed", color = "black") +
    geom_point(aes(color = Enrichment), alpha = 0.7, size = 2) +
    scale_color_manual(values = c(DOWN="darkcyan", NO="grey", UP="deeppink")) +
    geom_text_repel(
      data = label_data,
      aes(label = Genes, color = Enrichment),
      size = 4,
      max.overlaps = 50,
      fontface = "bold",
      show.legend = FALSE
    ) +
    labs(x = "log2 FC", y = "-log10 Q-value", title = title) +
    coord_cartesian(xlim = c(-3, 6), ylim = c(0, 8))
}

# 8) Filter to what you actually want to plot (massive speedup)
df_small <- df %>%
  filter(`Comparison (group1/group2)` %in% comparisons_keep,
         PTM.ModificationTitle %in% ptm_keep) %>%
  select(
    `Comparison (group1/group2)`,
    PTM.ModificationTitle,
    `AVG Log2 Ratio`,
    Qvalue,
    UniProtIds,
    Genes,
    ProteinGroups
  ) %>%
  prep_for_volcano()

# 9) Summary table (helps you decide what is worth interpreting)
summary_tbl <- df_small %>%
  group_by(`Comparison (group1/group2)`, PTM.ModificationTitle) %>%
  summarise(
    n = n(),
    n_up = sum(Enrichment == "UP"),
    n_down = sum(Enrichment == "DOWN"),
    n_sig = sum(Enrichment != "NO"),
    n_poi_sig = sum(Label == "poi" & Enrichment != "NO"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))

write.csv(summary_tbl, file.path(out_tab, "summary_counts_priority_PTMs.csv"), row.names = FALSE)

# 10) Generate plots
for (cmp in comparisons_keep) {
  for (ptm in ptm_keep) {
    
    d <- df_small %>%
      filter(`Comparison (group1/group2)` == cmp,
             PTM.ModificationTitle == ptm)
    
    if (nrow(d) < 200) next
    
    p <- make_volcano(d, paste(cmp, "|", ptm))
    
    cmp_safe <- stringr::str_replace_all(cmp, "[^A-Za-z0-9]+", "_")
    ptm_safe <- stringr::str_replace_all(ptm, "[^A-Za-z0-9]+", "_")
    
    outfile <- file.path(out_fig, sprintf("volcano_%s__%s__q%s_fc%s_POI.png",
                                          cmp_safe, ptm_safe, q_cut, log2fc_cut))
    
    ggsave(outfile, p, width = 32, height = 20, units = "cm", dpi = 300)
  }
}

message("Done. See: ", out_fig, " and ", out_tab)
