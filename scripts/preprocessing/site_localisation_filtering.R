# =============================================================
# Step 2: Quality filtering of PTM candidates
library(dplyr)
library(data.table)
library(ggplot2)

# --- 1. Load both files ---
ptm_candidates <- fread("data_raw/PTM1.tsv")
site_loc       <- fread("data_raw/ptm_open_search_site_localisation.tsv")

# --- 2. Build site whitelist from localisation report ---
# Summarise PTM.SiteProbability across all replicates
# and conditions for each unique site

site_whitelist <- site_loc %>%
  group_by(PG.UniProtIds, PTM.ModificationTitle) %>%
  summarise(
    N_observations   = n(),
    Mean_SiteProb    = mean(PTM.SiteProbability, na.rm = TRUE),
    Max_SiteProb     = max(PTM.SiteProbability, na.rm = TRUE),
    N_highconf_obs   = sum(PTM.SiteProbability >= 0.90, 
                           na.rm = TRUE),
    N_standard_obs   = sum(PTM.SiteProbability >= 0.75, 
                           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    MatchKey         = paste(PG.UniProtIds, 
                             PTM.ModificationTitle, sep = "_"),
    passes_standard  = Mean_SiteProb >= 0.75,
    passes_highconf  = Mean_SiteProb >= 0.90
  )

cat("Whitelist built:", nrow(site_whitelist), "unique sites\n")
cat("Passing standard (>=0.75):", 
    sum(site_whitelist$passes_standard), "\n")
cat("Passing high confidence (>=0.90):", 
    sum(site_whitelist$passes_highconf), "\n")

# --- 3. Prepare PTM candidates ---
# Extract UniProt ID from Group column and build MatchKey
# Remove contaminants (Cont_ prefix in Group)

n_raw <- nrow(ptm_candidates)

ptm_candidates <- ptm_candidates %>%
  mutate(
    UniProt_extracted = sub("_.*", "", Group),
    MatchKey = paste(UniProt_extracted, 
                     PTM.ModificationTitle, sep = "_")
  )

# --- 4. Apply filters in stages ---

# Stage 1: Remove contaminants
ptm_f1 <- ptm_candidates %>%
  filter(!grepl("^Cont", UniProt_extracted))
n_f1 <- nrow(ptm_f1)

# Stage 2: Evidence filters
ptm_f2 <- ptm_f1 %>%
  filter(`# Unique Total Peptides` >= 2,
         `# of Ratios`             >= 2)
n_f2 <- nrow(ptm_f2)

# Stage 3: Standard localisation filter (>=0.75)
ptm_standard <- ptm_f2 %>%
  filter(MatchKey %in%
           site_whitelist$MatchKey[site_whitelist$passes_standard])
n_standard <- nrow(ptm_standard)

# Stage 4: High confidence localisation filter (>=0.90)
ptm_highconf <- ptm_f2 %>%
  filter(MatchKey %in%
           site_whitelist$MatchKey[site_whitelist$passes_highconf])
n_highconf <- nrow(ptm_highconf)

# --- 5. Filtering summary ---
filter_summary <- data.frame(
  Stage = c(
    "Raw dataset",
    "Contaminants removed",
    "Evidence filtered\n(peptides>=2, ratios>=2)",
    "Localisation >=0.75",
    "Localisation >=0.90"
  ),
  N_rows = c(n_raw, n_f1, n_f2, n_standard, n_highconf)
) %>%
  mutate(
    N_removed    = c(0, 
                     n_raw    - n_f1,
                     n_f1     - n_f2,
                     n_f2     - n_standard,
                     n_standard - n_highconf),
    Pct_retained = round(N_rows / n_raw * 100, 1),
    Stage        = factor(Stage, levels = Stage)
  )

print(filter_summary)
write.csv(filter_summary,
          "results/tables/step2_filtering_summary.csv",
          row.names = FALSE)

# --- 6. Waterfall plot ---
p_waterfall <- ggplot(filter_summary,
                      aes(x = Stage, y = N_rows)) +
  geom_bar(stat = "identity", fill = "#2C7BB6",
           colour = "black", width = 0.6) +
  geom_text(aes(label = paste0(N_rows, "\n(",
                               Pct_retained, "%)")),
            vjust = -0.3, size = 3) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "PTM candidates: filtering stages",
       x = NULL,
       y = "Number of rows retained") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("results/figures/step2_filter_waterfall.png", p_waterfall,
       width = 8, height = 5, dpi = 300)

# --- 7. Modification survival summary ---
mod_survival <- ptm_standard %>%
  group_by(PTM.ModificationTitle) %>%
  summarise(
    N_rows            = n(),
    N_unique_proteins = n_distinct(Genes),
    N_unique_sites    = n_distinct(Group),
    N_comparisons     = n_distinct(`Comparison (group1/group2)`),
    .groups = "drop"
  ) %>%
  arrange(desc(N_unique_proteins))

print(mod_survival)
write.csv(mod_survival,
          "results/tables/step2_modification_survival_summary.csv",
          row.names = FALSE)

# --- 8. Export filtered datasets ---
write.csv(ptm_standard,
          "data/processed/PTM1_filtered_standard.csv",
          row.names = FALSE)

write.csv(ptm_highconf,
          "data/processed/PTM1_filtered_highconf.csv",
          row.names = FALSE)

message("Step 2 complete.")
message("Standard filtered dataset: ", n_standard, " rows")
message("High confidence dataset: ", n_highconf, " rows")

# =============================================================
#Artefact and fixed modification removal
# =============================================================

# Modifications to exclude - artefacts, fixed modifications,
# non-enzymatic, or non-biological in this context
artefact_ptms <- c(
  "Carbamidomethyl (C)",  # Fixed mod - sample prep alkylation
  "Oxidation (M)",         # Predominantly sample prep artefact
  "Didehydro",             # Dehydration artefact
  "Methyl (E)",            # Non-enzymatic, low biological relevance
  "Oxidation (P)",         # Sample prep artefact
  "Cation:Na (DE)",        # Sodium adduct - ionisation artefact
  "Cation:Mg[II]",         # Metal adduct - ionisation artefact
  "Dimethyl (KR)",         # Often artefactual
  "Formyl",                # Non-enzymatic artefact
  "Carbamylation (KR)",    # Urea artefact from sample prep
  "Methyl (KR)",           # Low confidence biological relevance
  "Hydroxymethyl",         # Non-enzymatic
  "Propionamide",          # Acrylamide adduct - sample prep
  "Carboxy",               # Non-enzymatic
  "Dioxidation (MW)",      # Oxidation artefact
  "Carboxymethyl",         # Iodoacetate adduct
  "Trioxidation (C)",      # Oxidation artefact
  "Cysteinyl",             # Non-enzymatic
  "Phosphoadenosine",      # Very low detection, unusual
  "Biotin"                 # Exogenous label not relevant here
)

# Biological PTMs to retain
biological_ptms <- c(
  "Phospho (STY)",
  "GlyGly (K)",
  "Acetyl (K)",
  "HexNAc (ST)",
  "Trimethyl (K)",
  "Crotonyl",
  "Myristoyl",
  "Hex (K)",
  "Sulfo (STY)",
  "Nitrosyl",
  "Oxidation (W)"  # Tryptophan oxidation has some biological
  # relevance unlike Met oxidation
)

# Filter to biological PTMs only
ptm_biological <- ptm_standard %>%
  filter(PTM.ModificationTitle %in% biological_ptms)

cat("Rows after biological PTM filter:", 
    nrow(ptm_biological), "\n")
cat("Unique proteins retained:", 
    n_distinct(ptm_biological$Genes), "\n")
cat("Unique sites retained:", 
    n_distinct(ptm_biological$Group), "\n")

# Summary by modification type after biological filter
mod_bio_summary <- ptm_biological %>%
  group_by(PTM.ModificationTitle) %>%
  summarise(
    N_rows            = n(),
    N_unique_proteins = n_distinct(Genes),
    N_unique_sites    = n_distinct(Group),
    .groups = "drop"
  ) %>%
  arrange(desc(N_unique_proteins))

print(mod_bio_summary)
write.csv(mod_bio_summary,
          "results/tables/step2_biological_ptm_summary.csv",
          row.names = FALSE)

write.csv(ptm_biological,
          "data/processed/PTM1_biological.csv",
          row.names = FALSE)

