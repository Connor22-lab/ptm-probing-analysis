rm(list=ls())

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
setwd("C:/Users/conno/Desktop/Final Year School/Final Year Project/data_raw")
outdir <- "C:/Users/yourname/Desktop/PTM_analysis/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

safe_comp <- gsub("[^A-Za-z0-9]+", "_", one_comp)
# ---- inputs ----
cand_path <- "ptm_closed_search_noimputation.tsv"   # your candidates export
one_comp  <- "Unstim / LPS_Nigericin"  # must match exactly

# Provide POIs as a character vector OR load from a file
poi_genes <- c("poi_gene_list")

# ---- load ----
cand <- fread(cand_path, sep="\t", header=TRUE)

# ---- filter to one comparison (comparison-wide; no POI filter) ----
df <- cand %>%
  filter(`Comparison (group1/group2)` == one_comp) %>%
  mutate(
    # Clean up gene field: take first symbol if multiple; adjust if your delimiter differs
    Gene1 = str_split(as.character(Genes), pattern="[,; ]+") |> sapply(\(x) x[1]),
    is_poi = Gene1 %in% poi_genes,
    minusLog10P = -log10(Pvalue),
    BH_q = p.adjust(Pvalue, method = "BH"),
    sig = BH_q < 0.05
  )

if(nrow(df) == 0){
  stop("No rows matched one_comp. Run sort(unique(cand$`Comparison (group1/group2)`)) and copy the exact string.")
}

# ---- choose what to label ----
# Label ALL POIs, plus optionally top non-POI hits to make plot informative
df <- df %>%
  mutate(
    label = ifelse(is_poi, Gene1, NA_character_)
  )
table(df$is_poi)

# ---- volcano plot ----
p <- ggplot(df, aes(x = `AVG Log2 Ratio`, y = minusLog10P)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_point(data = df %>% filter(sig), alpha = 0.8, size = 1.2) +
  geom_point(data = df %>% filter(is_poi), size = 1.5) +
  ggrepel::geom_text_repel(
    data = df %>% filter(is_poi),
    aes(label = label),
    max.overlaps = Inf,
    size = 3,
    box.padding = 0.3
  ) +
  theme_bw() +
  labs(
    title = paste("Volcano:", one_comp),
    subtitle = paste0("BH-FDR within this comparison; N = ", nrow(df),
                      " tests; min BH q = ", signif(min(df$BH_q, na.rm=TRUE), 3)),
    x = "AVG Log2 Ratio",
    y = "-log10(Pvalue)"
  )

print(p)
ggsave(
  filename = file.path(outdir, paste0("volcano_", safe_comp, ".png")),, p, width = 7.5, height = 5.5, dpi = 300)

# ---- optional: export the annotated table ----
fwrite(
  df, 
       file = file.path(outdir, paste0("comparisonwide_", safe_comp, "_with_BHq.tsv")),
  sep = "\t")

