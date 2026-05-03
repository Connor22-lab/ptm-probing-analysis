# PTM Probing Analysis

Final year project investigating post-translational modifications (PTMs) in the context of innate immunity using mass spectrometry open and closed search data.

## Overview

This project analyses PTM data from Spectronaut searches to identify differentially modified sites across conditions. The pipeline covers preprocessing, protein abundance correction, and statistical analysis with volcano plots.

## Repository Structure

```
├── scripts/
│   ├── preprocessing/         # Data cleaning and filtering
│   ├── analysis/              # Statistical analysis scripts
│   ├── visualization/         # Plotting scripts
│   └── ptm_exploration.ipynb  # Exploratory analysis notebook
├── data/
│   └── processed/             # Filtered/corrected datasets
├── results/                   # Output tables and figures
├── structures/                # AlphaFold structural models
└── docs/                      # Project documentation
```

## Pipeline

### Preprocessing
1. `contaminant_removal_ptm_selection.R` — Remove contaminants, select PTM sites
2. `fdr_recalculation.R` — Recalculate FDR after filtering
3. `site_localisation_filtering.R` — Filter by site localisation probability
4. `protein_abundance_correction.R` — Correct PTM intensities for protein abundance changes

### Analysis
- `condition_stratified_analysis.R` / `condition_stratified_intermediate.R` — Stratified differential analysis by condition
- `ptm_volcano_pipeline.R` / `ptm_volcano_protein_normalised.R` — Volcano plot generation
- `ptm_distribution.R` — PTM class and site distribution
- `open_search_filtering_waterfall.R` — Open search result filtering with waterfall plots
- `fc_direction_flip.R` — Fold change direction analysis
- `ptm_class_screen.py` / `ptm_qc_and_summary.py` — PTM classification and QC

## Requirements

- **R** (≥ 4.0) with packages: `tidyverse`, `ggplot2`, `ggrepel`, `data.table`
- **Python** (≥ 3.8) with packages: `pandas`, `numpy`
- Raw data from Spectronaut (not included — placed in `data_raw/`)
