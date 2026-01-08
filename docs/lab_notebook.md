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

## 2026-01-06 — Repo setup and logging baseline (LN-001)

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

2026-01-07 — PTM1 QC + stats diagnostics workflow (LN-002)

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