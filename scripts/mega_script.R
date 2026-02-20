# ============================================================
# SCRIPT: PTM volcano plots (protein-normalised, selective comps)
# BLOCK 1 — Setup: working dir, output dir, safe filenames
# ============================================================

rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)

# ----------------------------
# 1) Paths
# ----------------------------
wd <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw"

setwd(wd)

# 2) Run-level settings
# ----------------------------
IMPUTED <- FALSE   # TRUE = imputed run, FALSE = no imputation
# ----------------------------
# Output directory (depends on imputation)
# ----------------------------
base_outdir <- "C:/Users/conno/Desktop/Final Year School/Final Year Project/figures/volcano_protein_normalised"

outdir <- file.path(
  base_outdir,
  if (IMPUTED) "imputed" else "noImpute"
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# ----------------------------
# 3) Safe name helper (for files)
# ----------------------------
safe_name <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s*/\\s*", "_vs_", x)   # "A / B" -> "A_vs_B"
  x <- gsub("\\s+", "_", x)           # spaces -> underscores
  x <- gsub("[^A-Za-z0-9_\\-]+", "", x) # drop weird characters
  x
}

# ----------------------------
# 4) Plot filename helper
# ----------------------------
make_plot_filename <- function(comparison,
                               ptm_highlight = NULL,
                               normalised = TRUE,
                               imputed = IMPUTED,
                               ext = "png") {
  
  comp_tag <- safe_name(comparison)
  
  norm_tag <- if (isTRUE(normalised)) "protNorm" else "raw"
  imp_tag  <- if (isTRUE(imputed)) "imputed" else "noImpute"
  
  ptm_tag <- if (!is.null(ptm_highlight) && nzchar(ptm_highlight)) {
    paste0("_PTM-", safe_name(ptm_highlight))
  } else {
    ""
  }
  
  paste0("volcano_", comp_tag, "_", norm_tag, "_", imp_tag, ptm_tag, ".", ext)
}
# ============================================================
# BLOCK 2 — Load data + standardise key columns
# ============================================================

# ----------------------------
# 2.1) Input filenames (EDIT if needed)
# ----------------------------
file_prot <- if (IMPUTED) {
  "ptm_closed_search_PG_candidates.tsv"
} else {
  "ptm_closed_search_noimputation_PG_candidates.tsv"
}

file_ptm <- if (IMPUTED) {
  "ptm_closed_search.tsv"
} else {
  "ptm_closed_search_noimputation.tsv"
}

stopifnot(file.exists(file_prot))
stopifnot(file.exists(file_ptm))

# ----------------------------
# 2.2) Load tables
# ----------------------------
prot_raw <- fread(file_prot)
ptm_raw  <- fread(file_ptm)

# ----------------------------
# 2.3) Rename awkward / key columns immediately
# ----------------------------
prot <- prot_raw %>%
  rename(
    Comparison     = `Comparison (group1/group2)`,
    log2fc_prot    = `AVG Log2 Ratio`,
    p_prot         = Pvalue,
    q_prot         = Qvalue
  )

ptm <- ptm_raw %>%
  rename(
    Comparison     = `Comparison (group1/group2)`,
    log2fc_ptm     = `AVG Log2 Ratio`,
    p_ptm          = Pvalue,
    q_ptm          = Qvalue
  )

# ----------------------------
# 2.4) Normalise comparison string formatting
# ----------------------------
norm_comp <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s*/\\s*", " / ", x)
  x <- gsub("\\s+", " ", x)
  x
}

prot <- prot %>% mutate(Comparison = norm_comp(Comparison))
ptm  <- ptm  %>% mutate(Comparison = norm_comp(Comparison))

# ----------------------------
# 2.5) Quick sanity checks (non-fatal)
# ----------------------------
message("Protein comparisons: ", length(unique(prot$Comparison)))
message("PTM comparisons: ", length(unique(ptm$Comparison)))

# show any comparisons present in PTM but not protein
setdiff(unique(ptm$Comparison), unique(prot$Comparison))

# ============================================================
# BLOCK 3 — Define comparison subset (direction-agnostic)
# ============================================================

# Define pairs you want to analyse (directionless)
KEEP_PAIRS_MAIN <- list(
  c("Unstim", "P3C4"),
  c("Unstim", "LPS"),
  c("P3C4", "P3C4_Nigericin"),
  c("LPS",  "LPS_Nigericin")
)

INCLUDE_ENDPOINT    <- FALSE
INCLUDE_PRIMINGTYPE <- FALSE

KEEP_PAIRS_ENDPOINT <- list(
  c("Unstim", "P3C4_Nigericin"),
  c("Unstim", "LPS_Nigericin")
)

KEEP_PAIRS_PRIMINGTYPE <- list(
  c("P3C4", "LPS"),
  c("P3C4_Nigericin", "LPS_Nigericin")
)

KEEP_PAIRS <- c(
  KEEP_PAIRS_MAIN,
  if (INCLUDE_ENDPOINT) KEEP_PAIRS_ENDPOINT else list(),
  if (INCLUDE_PRIMINGTYPE) KEEP_PAIRS_PRIMINGTYPE else list()
)

# Build all acceptable "A / B" strings for matching your exports
pair_to_both_dirs <- function(a, b) c(paste(a, b, sep = " / "), paste(b, a, sep = " / "))

KEEP_COMPARISONS <- unique(unlist(lapply(KEEP_PAIRS, \(ab) pair_to_both_dirs(ab[1], ab[2]))))
KEEP_COMPARISONS <- norm_comp(KEEP_COMPARISONS)

# Check presence
present_ptm  <- sort(unique(ptm$Comparison))
present_prot <- sort(unique(prot$Comparison))

missing_in_ptm  <- setdiff(KEEP_COMPARISONS, present_ptm)
missing_in_prot <- setdiff(KEEP_COMPARISONS, present_prot)

# This now should be mostly empty; if not, it's a naming mismatch (typo/underscore/etc.)
if (length(missing_in_ptm) > 0) {
  message("NOTE: Some directions not present in PTM table (this is normal if only one direction exists):")
  print(missing_in_ptm)
}
if (length(missing_in_prot) > 0) {
  message("NOTE: Some directions not present in Protein table (normal if only one direction exists):")
  print(missing_in_prot)
}
# ============================================================
# BLOCK 4 — Filter to subset, flip to preferred direction,
#           merge protein FC, compute protein-normalised PTM FC
# ============================================================

# ----------------------------
# 4.1) Filter to just the comparisons of interest
# ----------------------------
ptm_sub  <- ptm  %>% filter(Comparison %in% KEEP_COMPARISONS)
prot_sub <- prot %>% filter(Comparison %in% KEEP_COMPARISONS)

#print what we’re actually keeping
message("PTM comparisons kept:")
print(sort(unique(ptm_sub$Comparison)))
message("Protein comparisons kept:")
print(sort(unique(prot_sub$Comparison)))
# ----------------------------
# 4.2) Define biological state ranking (for flipping direction)
# Unstim < Primed < Activated
# ----------------------------
rank_condition <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x == "Unstim" ~ 1,
    x %in% c("P3C4", "LPS") ~ 2,
    x %in% c("P3C4_Nigericin", "LPS_Nigericin") ~ 3,
    TRUE ~ NA_real_
  )
}

# ----------------------------
# 4.3) Flip to preferred direction using Condition Numerator/Denominator
# ----------------------------
flip_using_condition_cols <- function(df, fc_cols) {
  
  stopifnot(all(c("Condition Numerator", "Condition Denominator") %in% names(df)))
  
  df2 <- df %>%
    mutate(
      num = trimws(as.character(`Condition Numerator`)),
      den = trimws(as.character(`Condition Denominator`)),
      r_num = rank_condition(num),
      r_den = rank_condition(den),
      .flip = !is.na(r_num) & !is.na(r_den) & (r_num < r_den)
    )
  
  message("Rows flagged to flip: ", sum(df2$.flip), " / ", nrow(df2))
  
  # flip fold-changes where needed
  df2 <- df2 %>%
    mutate(across(all_of(fc_cols), ~ ifelse(.flip, -., .)))
  
  # swap the condition columns where needed
  df2 <- df2 %>%
    mutate(
      tmp = `Condition Numerator`,
      `Condition Numerator`   = ifelse(.flip, `Condition Denominator`, `Condition Numerator`),
      `Condition Denominator` = ifelse(.flip, tmp, `Condition Denominator`)
    ) %>%
    select(-tmp)
  
  # rebuild Comparison label from (possibly swapped) condition columns
  df2 <- df2 %>%
    mutate(Comparison = paste(`Condition Numerator`, `Condition Denominator`, sep = " / "))
  
  # drop helper cols
  df2 %>% select(-num, -den, -r_num, -r_den, -.flip)
}

# Apply to subsets
ptm_sub  <- flip_using_condition_cols(ptm_sub,  fc_cols = c("log2fc_ptm"))
prot_sub <- flip_using_condition_cols(prot_sub, fc_cols = c("log2fc_prot"))


# ----------------------------
# 4.4) Merge protein FC onto PTM rows + compute normalised FC
# Join key: Comparison + ProteinGroups
# ----------------------------
joined <- ptm_sub %>%
  left_join(
    prot_sub %>% select(Comparison, ProteinGroups, log2fc_prot, p_prot, q_prot),
    by = c("Comparison", "ProteinGroups")
  ) %>%
  mutate(
    log2fc_norm = log2fc_ptm - log2fc_prot
  )

# Optional: drop PTM rows where protein FC is missing (recommended for normalised volcanos)
joined <- joined %>% filter(!is.na(log2fc_prot))
# ---- ensure numeric columns are numeric (Spectronaut exports can be read as char) ----
joined <- joined %>%
  mutate(
    log2fc_ptm  = as.numeric(log2fc_ptm),
    log2fc_prot = as.numeric(log2fc_prot),
    log2fc_norm = as.numeric(log2fc_norm),
    p_ptm       = as.numeric(p_ptm),
    q_ptm       = as.numeric(q_ptm)
  )


# Final check: how many rows per comparison after merge?
merged_counts <- joined %>%
  count(Comparison, name = "n_rows") %>%
  arrange(desc(n_rows))

print(merged_counts)
# ============================================================
# BLOCK 5 — Plotting + exporting top hits per comparison
# - FILTERED volcanos: only PTM(s) of interest
# - GLOBAL volcanos: all PTMs, highlight PTM(s) of interest
# - Export Top N hits per comparison as CSV
# ============================================================

library(ggrepel)

# ----------------------------
# 5.1) Global settings (EDIT)
# ----------------------------
USE_SIGNIFICANCE <- "q"   # "p" or "q" (y-axis)
X_COL <- "log2fc_norm"    # "log2fc_norm" (recommended) or "log2fc_ptm"

Q_THRESH <- 0.05
P_THRESH <- 0.05
FC_THRESH <- 1

LABEL_TOP_N <- 15

# Toggle plot modes
MAKE_FILTERED_VOLCANOS <- TRUE   # PTM-only volcanos
MAKE_GLOBAL_VOLCANOS   <- TRUE   # all PTMs, highlight chosen PTM(s)

# PTMs of interest (exact match to PTM.ModificationTitle)
PTMS_OF_INTEREST <- c(
  # "Oxidation (M)",
   "Phospho (STY)",
   "GlyGly (K)",
   "Acetyl (K)"
)
# ----------------------------
# PTM colour map (EDITABLE)
# ----------------------------
PTM_COLOURS <- c(
  "GlyGly (K)"     = "forestgreen",
  "Acetyl (K)"     = "orange",
  "Phospho (STY)"  = "dodgerblue"
)

# fallback colour for PTMs not in the map
PTM_OTHER_COLOUR <- "grey70"

# POI highlighting (optional)
POI_UNIPROT <- c(
  # "P01584"
)
POI_GENES <- c(
  # "IL1B", "CYLD"
)

# In GLOBAL mode only: highlight POI across all PTMs, or only within PTM(s) of interest?
POI_ONLY_WITHIN_PTMS_OF_INTEREST <- TRUE

# Output
SAVE_PLOTS <- TRUE
PLOT_W <- 7
PLOT_H <- 6
PLOT_DPI <- 300

# Top hits CSV export
EXPORT_TOP_HITS <- TRUE
TOP_HITS_N1 <- 50
TOP_HITS_N2 <- 500

# ----------------------------
# 5.2) Helpers: significance + POI matching
# ----------------------------
get_sig_cols <- function(use = c("p", "q")) {
  use <- match.arg(use)
  if (use == "p") {
    list(val = "p_ptm", thr = P_THRESH, label = "-log10(p)")
  } else {
    list(val = "q_ptm", thr = Q_THRESH, label = "-log10(q)")
  }
}

# UniProtIds sometimes contain multiple IDs; this handles both single and delimited fields.
poi_matcher <- function(uniprot_col, gene_col, poi_uniprot, poi_genes) {
  hit_uni <- rep(FALSE, length(uniprot_col))
  hit_gene <- rep(FALSE, length(gene_col))
  
  if (!is.null(poi_uniprot) && length(poi_uniprot) > 0) {
    # match if any POI id appears as a whole token in UniProtIds (split on ; , space)
    pat <- paste0("(^|[; ,])(", paste(poi_uniprot, collapse = "|"), ")([; ,]|$)")
    hit_uni <- grepl(pat, uniprot_col)
  }
  
  if (!is.null(poi_genes) && length(poi_genes) > 0) {
    patg <- paste0("(^|[; ,])(", paste(poi_genes, collapse = "|"), ")([; ,]|$)")
    hit_gene <- grepl(patg, gene_col)
  }
  
  hit_uni | hit_gene
}

# ----------------------------
# 5.3) Helper: filename builders
# ----------------------------
make_plot_filename2 <- function(comparison, mode = c("filtered", "global"),
                                ptm_tag = NULL, x_col = X_COL, use_sig = USE_SIGNIFICANCE,
                                ext = "png") {
  mode <- match.arg(mode)
  comp_tag <- safe_name(comparison)
  norm_tag <- if (x_col == "log2fc_norm") "protNorm" else "raw"
  imp_tag  <- if (IMPUTED) "imputed" else "noImpute"
  sig_tag  <- paste0("y", use_sig)
  
  ptm_part <- if (!is.null(ptm_tag) && nzchar(ptm_tag)) paste0("_PTM-", safe_name(ptm_tag)) else ""
  paste0("volcano_", comp_tag, "_", mode, "_", norm_tag, "_", sig_tag, "_", imp_tag, ptm_part, ".", ext)
}

make_hits_filename <- function(comparison, n, ptm_tag = NULL, mode = c("filtered", "global")) {
  mode <- match.arg(mode)
  comp_tag <- safe_name(comparison)
  imp_tag  <- if (IMPUTED) "imputed" else "noImpute"
  ptm_part <- if (!is.null(ptm_tag) && nzchar(ptm_tag)) paste0("_PTM-", safe_name(ptm_tag)) else ""
  paste0("topHits_", comp_tag, "_", mode, "_", n, "_", imp_tag, ptm_part, ".csv")
}

# ----------------------------
# 5.4) Helper: volcano plot for a prepared df (already filtered if needed)
# mode:
#  - filtered: everything shown is PTM(s) of interest (class uses POI vs Significant vs Other)
#  - global: all PTMs shown; PTM(s) of interest highlighted; POI optional
# ----------------------------
plot_volcano_prepared <- function(d, comparison, mode = c("filtered", "global"),
                                  ptms_interest = PTMS_OF_INTEREST,
                                  x_col = X_COL, use_sig = USE_SIGNIFICANCE,
                                  fc_thresh = FC_THRESH, label_top_n = LABEL_TOP_N) {
  
  mode <- match.arg(mode)
  sig <- get_sig_cols(use_sig)
  
  # ---- compute stats ----
  d2 <- d %>%
    mutate(
     # clamp to avoid -log10(0) = Inf, and handle NA
sig_raw = .data[[sig$val]],
sig_clamped = pmax(sig_raw, 1e-300, na.rm = FALSE),
y = -log10(sig_clamped),
      pass_sig = .data[[sig$val]] <= sig$thr,
      pass_fc  = abs(.data[[x_col]]) >= fc_thresh,
      sig_hit  = pass_sig & pass_fc
    )
  
  # ---- PTM flags ----
  d2 <- d2 %>%
    mutate(
      is_ptm_interest = if ("PTM.ModificationTitle" %in% names(.)) {
        `PTM.ModificationTitle` %in% ptms_interest
      } else FALSE
    )
  
  # ---- PTM colour key ----
  d2 <- d2 %>%
    mutate(
      ptm_key = dplyr::case_when(
        mode == "global"   & is_ptm_interest ~ as.character(`PTM.ModificationTitle`),
        mode == "filtered"                  ~ as.character(`PTM.ModificationTitle`),
        TRUE                                ~ "Other"
      )
    )
  
  # ---- POI flag ----
  d2 <- d2 %>%
    mutate(
      poi_match = poi_matcher(
        uniprot_col = as.character(UniProtIds),
        gene_col    = as.character(Genes),
        poi_uniprot = POI_UNIPROT,
        poi_genes   = POI_GENES
      ),
      is_poi = if (mode == "global" && POI_ONLY_WITHIN_PTMS_OF_INTEREST && length(ptms_interest) > 0) {
        poi_match & is_ptm_interest
      } else {
        poi_match
      }
    )
  
  # ---- labels ----
  lab <- d2 %>%
    mutate(label = ifelse(!is.na(Genes) & Genes != "", Genes, ProteinGroups)) %>%
    mutate(label_me = is_poi) %>%
    bind_rows(
      d2 %>%
        arrange(desc(y)) %>%
        slice_head(n = label_top_n) %>%
        mutate(label = ifelse(!is.na(Genes) & Genes != "", Genes, ProteinGroups),
               label_me = TRUE)
    ) %>%
    distinct(Comparison, ProteinGroups, PTM.ModificationTitle, .keep_all = TRUE) %>%
    filter(label_me)
  
  # ---- build plot (NO dangling +) ----
  p <- ggplot(d2, aes(x = .data[[x_col]], y = y)) +
    geom_hline(yintercept = -log10(sig$thr), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
    
    # non-POI points
    geom_point(
      data = d2 %>% filter(!is_poi),
      aes(fill = ptm_key),
      shape = 21,
      size = 2.0,
      alpha = 0.85,
      colour = NA
    ) +
    
    # POI points (black outline)
    geom_point(
      data = d2 %>% filter(is_poi),
      aes(fill = ptm_key),
      shape = 24,
      size = 2.8,
      alpha = 0.95,
      colour = "black",
      stroke = 1
    ) +
    
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = label),
      size = 3,
      max.overlaps = Inf
    ) +
    
    scale_fill_manual(
      name = "PTM type",
      values = c(PTM_COLOURS, Other = PTM_OTHER_COLOUR),
      breaks = names(PTM_COLOURS)
    ) +
    
    labs(
      title = paste0("Volcano (", mode, "): ", comparison),
      x = if (x_col == "log2fc_norm")
        "Protein-normalised log2FC (PTM − Protein)"
      else
        "PTM log2FC",
      y = sig$label
    ) +
    
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# ----------------------------
# 5.5) Helper: export top hits
# Ranking: smallest p/q first (most significant), then largest |FC|
# ----------------------------
export_top_hits <- function(d, comparison, mode = c("filtered", "global"), n, ptm_tag = NULL,
                            x_col = X_COL, use_sig = USE_SIGNIFICANCE) {
  mode <- match.arg(mode)
  sig <- get_sig_cols(use_sig)
  
  out <- d %>%
    mutate(
      sig_val = .data[[sig$val]],
      abs_fc = abs(.data[[x_col]])
    ) %>%
    arrange(sig_val, dplyr::desc(abs_fc)) %>%
    slice_head(n = n)
  
  fname <- make_hits_filename(comparison, n = n, ptm_tag = ptm_tag, mode = mode)
  fwrite(out, file.path(outdir, fname))
}

# ----------------------------
# 5.6) Main loop: by comparison
# ----------------------------
comparisons_to_plot <- sort(unique(joined$Comparison))

for (comp in comparisons_to_plot) {
  
  df_comp <- joined %>% filter(Comparison == comp)
  
  # ===== GLOBAL volcano (all PTMs, highlight PTMs of interest) =====
  if (MAKE_GLOBAL_VOLCANOS) {
    p_global <- plot_volcano_prepared(
      d = df_comp,
      comparison = comp,
      mode = "global",
      ptms_interest = PTMS_OF_INTEREST
    )
    print(p_global)
    
    if (SAVE_PLOTS) {
      fname <- make_plot_filename2(comp, mode = "global",
                                   ptm_tag = if (length(PTMS_OF_INTEREST) == 1) PTMS_OF_INTEREST[1] else NULL)
      ggsave(file.path(outdir, fname), p_global, width = PLOT_W, height = PLOT_H, dpi = PLOT_DPI)
    }
    
    if (EXPORT_TOP_HITS) {
      export_top_hits(df_comp, comp, mode = "global", n = TOP_HITS_N1,
                      ptm_tag = if (length(PTMS_OF_INTEREST) == 1) PTMS_OF_INTEREST[1] else NULL)
      export_top_hits(df_comp, comp, mode = "global", n = TOP_HITS_N2,
                      ptm_tag = if (length(PTMS_OF_INTEREST) == 1) PTMS_OF_INTEREST[1] else NULL)
    }
  }
  
  # ===== FILTERED volcanos (only PTM(s) of interest) =====
  if (MAKE_FILTERED_VOLCANOS) {
    
    # If no PTMs specified, we skip (filtered mode requires PTM list)
    if (length(PTMS_OF_INTEREST) == 0) {
      message("Filtered mode skipped for ", comp, " because PTMS_OF_INTEREST is empty.")
      next
    }
    
    # Option 1: one filtered volcano per PTM (recommended for clarity)
    for (ptm_name in PTMS_OF_INTEREST) {
      
      df_filt <- df_comp %>% filter(PTM.ModificationTitle == ptm_name)
      
      # If a comparison has zero rows for that PTM, skip gracefully
      if (nrow(df_filt) == 0) {
        message("No rows for ", ptm_name, " in ", comp, " (filtered volcano skipped).")
        next
      }
      
      p_filt <- plot_volcano_prepared(
        d = df_filt,
        comparison = comp,
        mode = "filtered",
        ptms_interest = ptm_name  # not used much in filtered mode, but fine
      )
      print(p_filt)
      
      if (SAVE_PLOTS) {
        fname <- make_plot_filename2(comp, mode = "filtered", ptm_tag = ptm_name)
        ggsave(file.path(outdir, fname), p_filt, width = PLOT_W, height = PLOT_H, dpi = PLOT_DPI)
      }
      
      if (EXPORT_TOP_HITS) {
        export_top_hits(df_filt, comp, mode = "filtered", n = TOP_HITS_N1, ptm_tag = ptm_name)
        export_top_hits(df_filt, comp, mode = "filtered", n = TOP_HITS_N2, ptm_tag = ptm_name)
      }
    }
  }
}
# ============================================================
# BLOCK 6 — Save flipped + normalised datasets
# ============================================================

SAVE_DATASETS <- TRUE

if (SAVE_DATASETS) {
  
  tag_imp <- if (IMPUTED) "imputed" else "noImpute"
  
  # PTM-only (flipped)
  fwrite(ptm_sub,  file.path(outdir, paste0("ptm_flipped_", tag_imp, ".csv")))
  
  # Protein-only (flipped)
  fwrite(prot_sub, file.path(outdir, paste0("protein_flipped_", tag_imp, ".csv")))
  
  # Joined (flipped + protein-normalised)
  fwrite(joined,   file.path(outdir, paste0("ptm_joined_protNorm_", tag_imp, ".csv")))
  
  message("Saved datasets into: ", outdir)
}

