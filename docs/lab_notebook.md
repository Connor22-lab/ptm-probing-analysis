# Laboratory Notebook — Final Year Project

**Project title:** Computational discovery on novel post-translational modifications in the inflammasome system  
**Student:** Connor  
**Notebook format:** Markdown (version-controlled, reproducible)  
**Primary data type:** Spectronaut PTM probing exports (e.g., `PTM1.tsv`)  
**Last updated:** 2026-01-06 (Europe/Dublin)

---

## How this notebook is organised

- **Chronological entries** (LN-###): one entry per work session or meeting.
- Each entry records:
  - objective → inputs → actions → outputs → decisions → issues → next steps.
- Code is **referenced** (scripts and filenames) rather than pasted in full.
- “Intermediate outputs go to results/ (ignored by Git). Curated deliverables go to figures/ and tables/ (tracked).”
## Quick links
- [Log entries](#log-entries)

---

## Ethics and appropriate use of AI

### Ethics statement (computational reanalysis)
This project involves computational reanalysis of existing, anonymized proteomics data from the THP-1 immortalized cell line. No new biological samples are collected. Primary considerations are data management and research integrity, addressed by:
- Storing/processing data on secure, university-approved systems.
- Fully documenting methods for reproducibility.

### Appropriate use of AI (agreed)
AI may be used for:
- Identifying research gaps and supporting literature review.
- Methodology planning and design.
- Debugging code.
- Finding relevant papers for discussion support.

AI must **not** be used for dissertation text generation.

---

## Project structure (recommended)

Final Year Project/
  data_raw/            # Spectronaut exports (ignored by Git)
  data_processed/      # cleaned intermediates (ignored by Git)
  scripts/             # python + r scripts (numbered/runnable)
  results/             # scratch outputs (ignored by Git)
  figures/             # curated figures (tracked)
  tables/              # curated tables (tracked)
  docs/
    lab_notebook.md
    archive/
      retrospective_reconstruction.md
    decisions.md
    meeting_notes.md

```

---

## Retrospective summary (pre-notebook work)

This summary records key background work completed before converting to a formal notebook format.

### High-level aim
Discover and prioritise novel PTM sites in inflammasome signalling, using open-search/PTM discovery and Spectronaut PTM probing, then contextualise prioritised sites structurally.

### Timeline (key milestones)
- **2025-09-29:** First supervisor meeting; project introduced; reading and onboarding resources provided.
- **2025-10-03:** Follow-up supervisor meeting; limited progress due to illness.
- **2025-10-23:** Drafted rough research proposal.
- **2025-10-30:** Submitted research proposal.
- **2025-11-19:** Supervisor email: Spectronaut expected as primary analysis platform; advised Spectronaut tutorials/playlists.
- **2025-11-27:** Introduced to postdoc for Spectronaut/PTM workflow support.
- **2025-12-03:** Met postdoc; walked through Spectronaut PTM workflow; initiated first run; generated initial export (`PTM1.tsv` / PTM probing output).
- **2025-12-16:** Reviewed PTM workflow outputs; agreed to prioritise volcano plots and candidate lists; plan to track PTMs across inflammasome-relevant conditions.

---

## Log entries

## LN-### — Template for future entries

**Date / Time:**  
**Objective:**  
**Inputs (files):**  
**Environment:**  
**Methods / Actions:**  
**Results / Outputs (filenames):**  
**Interpretation / Decisions:**  
**Issues / Fixes:**  
**Next steps:**  
**Meeting notes (if applicable):**  
**AI use (if any; allowed uses only):**  

## LN-001 
2026-01-06 — Repo setup and logging baseline 

**Objective:** Initialise Git/GitHub repo and establish reproducible documentation + file governance.

**Inputs (files):** N/A (repo setup).

**Environment:** VS Code on Windows; Git configured with university email.

**Methods / Actions:**
- Initialised Git repository in `Final Year Project/`.
- Added `.gitignore` to exclude `data_raw/`, `results/`, `data_processed/`, and Office temp files.
- Created `docs/archive/retrospective_reconstruction.md` containing reconstructed LN-001–LN-009.
- Updated `docs/lab_notebook.md` to include compressed background + live entries.

**Results / Outputs (filenames):**
- `.gitignore`
- `docs/lab_notebook.md`
- `docs/archive/retrospective_reconstruction.md`

**Interpretation / Decisions:**
- Deliverables saved to `figures/` and `tables/`; intermediates saved to `results/`.

**Issues / Fixes:** None.

**Next steps:**
- Load `PTM1.tsv` from `data_raw/` in Python and generate QC summary (shape, missingness, Qvalue distribution, #Ratios distribution).
- Generate first PTM summary table per Comparison × PTM type; export curated summary to `tables/`.

**AI use:** Used AI to create and format repo

## LN-002
2026-01-07 — PTM1 QC + stats diagnostics workflow 

**Objective:** Set up a reproducible Python workflow to load data_raw/PTM1.tsv, generate QC summaries/plots, and export a curated PTM summary table (Comparison × PTM type). Validate whether Pvalue/Qvalue behave as expected for downstream inference.

**Inputs (files):** data_raw/PTM1.tsv

**Environment:** Windows + VS Code, Conda env: ptm-qc (Python) Key packages: pandas, numpy, matplotlib

**Methods / Actions:**
- Configured VS Code to use the ptm-qc interpreter and verified the active Python executable.
- Implemented a reproducible analysis script to:
  load PTM1.tsv into pandas,
  compute dataset shape and column-level missingness,
  plot distributions for Qvalue, -log10(Qvalue), and # of Ratios,
  aggregate to a summary table by Comparison (group1/group2) × PTM.ModificationTitle, including counts of significant hits (by Qvalue threshold).
  Implemented a diagnostics script (qvalue_check.py) to:
  list p/q-related columns,
  quantify distributions of Pvalue and Qvalue,
  compute overall and per-comparison counts passing common Qvalue thresholds,
  compare nominal p-value signal versus FDR-adjusted discoveries.

**Results / Outputs (filenames):**
 - tables/PTM1_missingness_by_column.csv
 - tables/PTM1_summary_by_comparison_ptmtype.csv
 - figures/PTM1_qvalue_distribution.png
 - figures/PTM1_minuslog10_qvalue_distribution.png
 - figures/PTM1_num_ratios_distribution.png

**Script(s) created/updated:**
 - scripts/ptm_qc_and_summary.py
 - scripts/qvalue_check.py

**Interpretation / Decisions:**
 - Quantification completeness appears reasonable (most rows have ratios; minimal missingness in core columns).
 - There is a substantial nominal p-value tail (many rows with P ≤ 0.05), but almost no rows pass Q ≤ 0.05 (only single-digit discoveries across ~828k rows).
 - Interpreted as a strong multiple-testing burden and/or overly granular export level; FDR becomes extremely stringent at this scale.

Decision: treat this run as a QC/methods checkpoint; defer biological interpretation/volcano plots until the Spectronaut export/statistical level is revised (e.g., site-level aggregation and/or filtering).

**Issues / Fixes:**
Confirmed VS Code interpreter/env selection to ensure correct package availability and reproducible runs.

**Next steps:**
 - When Spectronaut access is available: review analysis configuration (replicate/group assignment, statistical testing settings) and export level (feature vs site vs protein) to reduce the number of tests and obtain interpretable FDR-controlled results.
 - Update summary reporting to include effect-size and nominal p-value counts (e.g., P ≤ 0.05 and |log2FC| thresholds) for triage and comparison-level prioritisation prior to re-export.

**Meeting notes (if applicable):**
 - None.

**AI use (if any; allowed uses only):**
Used AI assistance to (i) draft Python scripts for QC and summarisation, (ii) troubleshoot Windows/VS Code/conda execution issues, and (iii) interpret diagnostic outputs (Pvalue vs Qvalue behaviour under multiple testing).

## LN-003 — 

2026/01/10
**Objective:** 
Prioritise biologically relevant PTM class for follow-up Spectronaut Analysis

**Inputs (files):**  PTM1.tsv
scripts/ptm_class_screen.py (created in VS Code)

**Environment:**  VS Code, python

**Methods / Actions:** Implemented a PTM class screen in python to identify PTM classes
Evidence Filtering:
# of Ratios >= 2

# Unique Total Peptides >= 2

Absolute AVG Log2 Ratio >= 0.3

**Results / Outputs (filenames):**  
results_small/ptm_class_summary.csv
results_small/ptm_class_summary_by_comparison.csv

**Interpretation / Decisions:**  
PTM_class_summary.csv showed candidates: Phospho, GlyGly, HexNAc, Acetyl

**Issues / Fixes:**  
N/A

**Next steps:** 
Use ptm_class_summary_by_comparison.csv to confirm which PTM class is strongest in the primary biological comparison (e.g., priming vs activation contrast).

Extend the script to automatically export Top N Phospho candidates per comparison (rank by Pvalue, then absolute effect size, then evidence/support).

In Spectronaut, perform a focused re-analysis for Phospho (STY) (avoid wide PTM probing if possible) and export site-resolved PTM outputs (PTM site report) plus an optional run-pivot table to assess missingness/replicate consistency. 

**Meeting notes (if applicable):**  
N/A
**AI use (if any; allowed uses only):** Used AI to create ptm_class_screen.py

## LN-004 — 

**Date:** 2026-01-13

**Objective**
Generate interpretable volcano plots for inflammasome state transitions, focusing on Activation (Active vs Primed) for two priming conditions (LPS and P3C4).

Validate whether the broad PTM probing export yields statistically meaningful candidates.

Decide next analytical step and whether a narrower PTM search in Spectronaut is required to increase power.

**Inputs (files):** PTM1.tsv and POI.csv  
**Environment:** RStudio 
**Methods / Actions:** 
Comparison selection
selected activation contrasts for each priming condition:
LPS / LPS+Nigericin
P3C4 / P3C4+Nigericin

Created subsets:
Cands1_LPS_act
Cands1_P3C4_act

To make the x-axis biologically consistent (positive = higher in active state), implemented:
log2FC_state = -AVG Log2 Ratio
Plotted x-axis as log2FC_state and y-axis as -log10(Qvalue).

Volcano plot generation
Applied conventional thresholds:
log2FC ≥ 0.58 (≈1.5×), Q < 0.05
Annotated points as UP/DOWN/NO based on these thresholds.
Generated unlabelled, labelled, and POI-labelled volcano plots for both LPS and P3C4 activation.

POI mapping
Loaded POI UniProt list and implemented pattern-based matching for UniProtIds (to account for multiple IDs per row) so POI labelling was not dependent on exact single-ID matches. 

**Results / Outputs (filenames):** 
Candidate tables:
Candidates_LPS_act_tr1.csv
Candidates_P3C4_act_tr1.csv

Volcano plots:
VolcanoPlot1.png through VolcanoPlot10.png (variants with/without labels and POIs) 

Volcano plots for both activation contrasts show no features meeting Q < 0.05 with the chosen effect-size threshold, consistent with the earlier finding that the PTM probing was too broad to produce statistically robust hits after multiple testing correction.

The point clouds are heavily concentrated at low -log10(Qvalue) with a symmetric spread around log2FC ≈ 0, indicating limited evidence for strong differential PTM signals at the current search breadth.

POI genes (e.g., CASP1, PYCARD, IL1B, GSDMD) appear in the datasets and can be labelled, but they remain low-significance under the current global multiple-testing burden.

**Interpretation:** 
The absence of significant hits is likely driven by a combination of:
Large hypothesis space (many PTM classes / features tested),
Limited power per feature (variable peptide evidence / replicate counts),
Multiple-testing correction penalising broad probing.

This supports the conclusion that the current “wide PTM probing” approach is better treated as exploratory and requires narrowing to generate actionable candidate PTMs. 

**Decision:** Rerun the analysis in Spectronaut with a narrower PTM search.

Rationale: Reducing the number of PTM classes and the overall feature space should:
Increase statistical power
Reduce the multiple-testing burden (improving Q-values),
Yield clearer biological interpretation aligned with inflammasome activation biology.

**Issues / Fixes:** 
Minor changes to script in relation to labeling and x and y axis values
 
**Next steps:** 
Restrict PTMs to a small, biologically motivated set (example starting set):
Phospho (STY) (signalling/kinase activity)
Acetyl (K) (regulatory)
GlyGly (K) (ubiquitination proxy)
HexNAc (ST) (O-GlcNAc / glyco-related signalling)
Maintain the same comparison structure (Activation contrasts) for direct comparability with the current outputs.
 
**Meeting notes (if applicable):** N/A 

**AI use (if any; allowed uses only):** 
Used to assist changes in script, relating to log2FCstate

## LN-005 --

**Date:**  16/01/2026

**Objective:** Generate Volcano Plots for different Comparisons and Create pivot tables for genes of interest

**Environment:**  R, Excel

**Methods / Actions:** 
Merge PTM1.tsv with genes of interest list, create pivot table from merged dataset, add slicer to filter by PTM modification
Ran Rscripts to generate a series of volcano plots, done comparisons for LPS and P3C4 across different conditions

**Results / Outputs (filenames):**  Series of volcano plots, ptms_per_gene,xlsx and revised_ptm_per_gene.xlsx

**Interpretation / Decisions:** Need to run a closed search as very few modifications of signfigance

**Issues / Fixes:** created a more robust gene list as was missing potentially important proteins with modifications

**Next steps:** Present preliminary results and run a closed search 

**Meeting notes (if applicable):**  N/A

**AI use (if any; allowed uses only):** help with excel merge formatting

## LN-006

**Date:** 18/01/2026 
**Objective:**  Generate a series of csv and volcano plots for different comparisons to dictate closed search, noting POI and adding them to pivot table/notes
**Inputs (files):**  PTM1.tsv , ptm_highlight_volcano.R
**Environment:**  R
**Methods / Actions:** Ran script changing comparison and modification type, added line to display all comparisons to reduce error  
**Results / Outputs (filenames):**  too many to list
**Interpretation / Decisions:**  Hard to interpret an open search but some intersting proteins were modified that are linked to inflammasome regulation
**Issues / Fixes:** N/A  
**Next steps:** Communicate results with DINO and decide next step
**Meeting notes (if applicable):**  N/A
**AI use (if any; allowed uses only):** none

## LN-007

**Date:** 23/01/26   

**Objective:**  Investigate Data from pivot table, create flow chart

**Inputs (files):**  PTM1.tsv, revised_ptm_per_gene.xlsx
**Environment:**  Excel, RStudio, PowerPoint

**Methods / Actions:**  The values for PTMs detected on select proteins were very unusual.
After selecting a few from the table, such as IL1B and HSP, I decided to cross check all the columns.
Discovered that proteins with 15 ptms,there was one per comparison, all at the same site.
Follow up diagnostics in R to determine similarity in modifications across proteins.
Workflow creation in Powerpoint, detailing the MS DATA, Spectronaut search, R work, observations and next steps


**Results / Outputs (filenames):**  
PTMs were nearly always consistent across conditions, especially with inflammasome components.
PTM Workflow.pptx

**Interpretation / Decisions:**  
Data has some serious flaws in strength, good for a prelimnary search and to get used to the software and how to analyse the data sets.

**Issues / Fixes:** N/A 

**Next steps:** Start to write up my Visual Abstract, and presentation for bens, try to combine them 

**Meeting notes (if applicable):**  N/A

**AI use (if any; allowed uses only):** N/A

## LN-007 

**Date:** 2026/01/27

**Objective:** Explore Parameters on Spectronaut, Generate reports for site localisation, Test AlphaMap, evaluate closed search data

**Inputs (files):** ptm_closed_search.tsv, scripts/Script1_Bar graphs_all_ptms.R, scripts/Script2_volcano_ptm_sample.R

**Environment:**  Spectrnaut, AlphaMap, Excel

**Methods / Actions:**  
First I exported the list of candidates with and without global imputation.
Tested the effect of cross run normalisation.
Generated Report for AlphaMap
Generated Report for site localisation (adjusted the threshold to 0.5)
Imported site localisation tsv into excel, merged with POI table for ease of use
Imported closed search into R, ran script for bar charts and made slight changes Script2 to change OUTDIR and labelling of files.

**Results / Outputs (filenames):**results/tables/ptm_localisation_poi.xlsx, data_raw/ptm_closed_search_alphamap_export.tsv, data_raw/ptm_closed_search_noimputation.tsv, data_raw/ptm_closed_search_site_localisation.tsv, 

**Interpretation / Decisions:** Unsure of specificity/reliabiliy of site localisation in Spectronaut, lots of whole values/suspicous results.
not much of a difference between closed and open search, Q-values are all non-significant
May have to correct the FDR by removing fixed modification from statistical testing

**Issues / Fixes:** N/A
**Next steps:**  statistical testing, volcano plot generation of POI. evaluate site localisation data. work on visual abstract and presentation for ben
**Meeting notes (if applicable):**  N/A
**AI use (if any; allowed uses only):** troubleshooting in Excel.