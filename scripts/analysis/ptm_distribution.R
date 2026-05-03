rm(list=ls())

library(here)
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)

outdir <- here("results", "figures", "closed_search")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

site <- fread(here("data_raw", "ptm_closed_search_site_localisation.tsv"), sep="\t", header=TRUE)

# --- localisation filter (0.5 is permissive; 0.75 is more standard) ---
site <- site %>%
  mutate(PTM.SiteProbability = as.numeric(PTM.SiteProbability)) %>%
  filter(!is.na(PTM.SiteProbability), PTM.SiteProbability >= 0.5)

# --- deduplicate to UNIQUE sites per condition x PTM type ---
site_unique <- site %>%
  distinct(R.Condition, PTM.ModificationTitle, PTM.CollapseKey)

ptm_counts <- site_unique %>%
  count(R.Condition, PTM.ModificationTitle, name = "n_unique_sites")

# Optional: order conditions in a sensible biological order
ptm_counts <- ptm_counts %>%
  mutate(R.Condition = factor(R.Condition,
                              levels = c("Unstim","LPS","P3C4","Nigericin","LPS_Nigericin","P3C4_Nigericin")))

# ---------- helpers ----------
# Totals per condition (for labeling bar tops)
totals <- ptm_counts %>%
  group_by(R.Condition) %>%
  summarise(total_unique_sites = sum(n_unique_sites), .groups="drop")

ptm_palette <- c(
  "Phospho (STY)"        = "#1f77b4",  # blue
  "Acetyl (K)"           = "#ff7f0e",  # orange
  "GlyGly (K)"           = "#2ca02c",  # green
  "Oxidation (M)"        = "#D3D3D3",  # light grey
  "Carbamidomethyl (C)" = "#7f7f7f"   # grey
)
ptm_order <- c(
  "GlyGly (K)",
  "Acetyl (K)",
  "Phospho (STY)",
  "Oxidation (M)",
  "Carbamidomethyl (C)")

ptm_counts$PTM.ModificationTitle <- factor(
  ptm_counts$PTM.ModificationTitle,
  levels = ptm_order)

# =========================
# PLOT A: QC (ALL PTMs) - show totals only
# =========================
ymax <- max(totals$total_unique_sites)

p_count_all <- ggplot(ptm_counts,
                      aes(x = R.Condition, y = n_unique_sites, fill = PTM.ModificationTitle)) +
  geom_bar(stat="identity", width=0.8) +
  geom_text(
    data = totals,
    aes(x = R.Condition, y = total_unique_sites, label = comma(total_unique_sites)),
    vjust = -0.4,
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_fill_manual(values = ptm_palette, drop = FALSE) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Unique PTM sites",
    title = "Unique PTM sites per condition"
  ) +
  coord_cartesian(clip = "off")

print(p_count_all)
ggsave(file.path(outdir, "PTM_unique_sites_counts_ALL_QC.png"),
       plot = p_count_all, width = 8.5, height = 6.3, dpi = 300)

# =========================
# PLOT B: Biology view (DROP Carbamidomethyl + Oxidation)
# =========================
drop_ptms <- c("Carbamidomethyl (C)", "Oxidation (M)")

ptm_counts_bio <- ptm_counts %>%
  filter(!PTM.ModificationTitle %in% drop_ptms)

# Add proportions within condition (for percent stacked option)
ptm_props_bio <- ptm_counts_bio %>%
  group_by(R.Condition) %>%
  mutate(prop = n_unique_sites / sum(n_unique_sites)) %>%
  ungroup()

# Label segments only if they are big enough to avoid clutter
ptm_props_bio <- ptm_props_bio %>%
  mutate(label_n = ifelse(prop >= 0.05, comma(n_unique_sites), NA))

p_prop_bio <- ggplot(ptm_props_bio,
                     aes(x = R.Condition, y = prop, fill = PTM.ModificationTitle)) +
  geom_bar(stat="identity", width=0.8) +
  geom_text(
    aes(label = label_n),
    position = position_stack(vjust = 0.5),
    size = 3,
    na.rm = TRUE
  ) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = ptm_palette, drop = TRUE) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Proportion of unique PTM sites",
    title = "PTM type distribution across conditions",
    subtitle = "Carbamidomethyl/Oxidation removed"
  )

print(p_prop_bio)
ggsave(file.path(outdir, "PTM_unique_sites_stackedbar_BIO_proportion.png"),
       plot = p_prop_bio, width = 8.5, height = 5.5, dpi = 300)
