# Laboratory Notebook — Final Year Project

**Project title:** Computational discovery on novel post-translational modifications in the inflammasome system  
**Student:** Connor  
**Notebook format:** Markdown (version-controlled, reproducible)  
**Primary data type:** Spectronaut PTM probing exports (e.g., `PTM1.tsv`)  
**Last updated:** 2026-01-04 (Europe/Dublin)

---

## How this notebook is organised

- **Chronological entries** (LN-###): one entry per work session or meeting.
- Each entry records:
  - objective → inputs → actions → outputs → decisions → issues → next steps.
- Code is **referenced** (scripts and filenames) rather than pasted in full.
- Figures/tables are stored in `results/` and linked by filename.

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

Create and maintain the following structure in your project directory:

```
Final Year Project/
  data_raw/            # Spectronaut exports (e.g., PTM1.tsv)
  data_processed/      # filtered/cleaned subsets you generate
  scripts/             # R scripts (numbered, runnable in order)
  results/
    figures/
    tables/
  docs/
    lab_notebook.md    # this file
    decisions.md       # optional (global decisions)
    meetings.md        # optional (standalone meeting minutes)
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

## Global decisions (living list)

1) **Primary analytical framing:** condition-based contrasts that reflect inflammasome biology (e.g., priming vs activation) rather than timepoints (unless timepoints are present in a future dataset).
2) **Initial significance thresholds (for exploratory volcano plots):**
   - log2FC cutoff: ±0.58 (≈1.5-fold)
   - Q-value cutoff: 0.05  
   These may be tightened later depending on PTM class, localisation confidence, and replication.
3) **POI labelling strategy:** proteins of interest stored in `POI.csv` as UniProt accessions; matching must handle semicolon-delimited `UniProtIds` values from Spectronaut.
4) **Reproducibility:** all plots saved with informative filenames including comparison + PTM type + thresholds.

---

# Entries

## LN-001 — 2025-09-29 — Supervisor meeting (project onboarding)

**Objective:** Project introduction, define initial reading list and onboarding tasks.  
**Actions / notes:**
- Received project overview: inflammasome biology and role of PTMs as regulatory switches.
- Provided background reading and onboarding resources (review + NLRP3 priming/PTM paper; Spectronaut resources).
**Decisions:**
- Begin with literature familiarisation, then proceed to software onboarding.
**Next steps:**
- Read assigned core review(s).
- Begin basic familiarisation with Spectronaut and PTM probing concepts.

---

## LN-002 — 2025-10-03 — Supervisor follow-up meeting

**Objective:** Progress check-in.  
**Notes:**
- Limited progress due to illness.
**Next steps:**
- Resume reading and begin drafting research proposal outline.

---

## LN-003 — 2025-10-23 — Drafted rough research proposal

**Objective:** Create initial written framing (context, gap, aims, methods).  
**Outputs:**
- Draft proposal (working document): `Rough Research Plan_BC.docx` (local).
**Decisions:**
- Frame inflammasome responses across phases (priming, activation/assembly, termination) and treat PTMs as dynamic regulators.
**Next steps:**
- Finalise and submit proposal; align with module requirements (ethics + AI statement + progress plan).

---

## LN-004 — 2025-10-30 — Research proposal submitted

**Objective:** Submit project research proposal.  
**Outputs:**
- Submitted proposal PDF (local copy).
**Next steps:**
- Begin software onboarding and identify initial dataset(s) for PTM probing.

---

## LN-005 — 2025-11-19 — Supervisor guidance (Spectronaut focus)

**Objective:** Confirm primary analysis platform and onboarding pathway.  
**Notes:**
- Spectronaut likely primary platform; access via remote licensed computers.
- Suggested tutorial playlists for baseline familiarity.
**Resources recorded:**
- https://youtube.com/playlist?list=PLci6bfaTzAjRxryp8cABtNjYvGK0eeRBw&si=VShm5UVS0Ix4o5DW
- https://youtube.com/playlist?list=PLci6bfaTzAjQCguHZO6YBghn77TmluySh&si=ZNuIUkCuU3SpZlQv
**Next steps:**
- Watch core Spectronaut tutorials; prepare questions for postdoc session.

---

## LN-006 — 2025-11-27 — Introduced to postdoc for PTM workflow support

**Objective:** Establish technical support for Spectronaut PTM workflow.  
**Next steps:**
- Arrange onboarding session; prepare dataset and questions.

---

## LN-007 — 2025-12-03 — Postdoc onboarding (Spectronaut PTM workflow)

**Time:** 10:00–10:30  
**Objective:** Walk through Spectronaut PTM probing workflow; initiate first run.  
**Inputs:**
- Spectronaut software access (remote licensed environment).
- Dataset configured for PTM probing run.
**Actions:**
- Walked through PTM probing workflow.
- Started first run and exported initial output (recorded as `PTM1.tsv`).
**Outputs:**
- `PTM1.tsv` (Spectronaut export; stored under `data_raw/` once localised).
**Next steps:**
- Review export structure; determine comparison design; plan first volcano plots.

---

## LN-008 — 2025-12-16 — Review PTM workflow outputs and plotting plan

**Objective:** Decide initial analysis outputs suitable for early results.  
**Notes / actions:**
- Reviewed analysis plots within PTM workflow.
- Focused on volcano plots and selecting proteins/modifications of interest.
- Discussed exporting candidate lists (with/without global imputation) for custom plotting in R.
**Decisions:**
- First tangible results: modification-summary trends and volcano plots.
- Plan to generate candidate lists per comparison/PTM type and label inflammasome proteins of interest.
**Next steps:**
- Produce R plotting workflow; create POI list; choose primary comparisons aligned to priming/activation.

---

## LN-009 — 2026-01-04 — PTM1.tsv triage in R (data understanding)

**Objective:** Validate the structure of `PTM1.tsv` and identify the correct way to subset (no timepoints; condition comparisons).  
**Inputs:**
- `PTM1.tsv` (Spectronaut export; unmodified).
**Environment:**
- R with `data.table` (`fread`), plus plotting stack planned (`ggplot2`, `ggrepel`, `dplyr`, `stringr`).
**Actions:**
- Loaded `PTM1.tsv` and confirmed schema.
- Verified: 828,005 rows × 21 columns.
- Enumerated comparisons: 15 pairwise condition contrasts (Unstim, P3C4, LPS, Nigericin, and combinations).
- Confirmed required volcano columns exist: `AVG Log2 Ratio`, `Qvalue`, `Genes`, `UniProtIds`, and `PTM.ModificationTitle`.
**Key observation:**
- Dataset is **condition-based**, not time-course; therefore analysis should be separated by:
  1) `Comparison (group1/group2)` and
  2) `PTM.ModificationTitle`.
**Decisions:**
- Proceed by selecting a small set of biologically meaningful contrasts (priming vs activation framing) and generating volcano plots per PTM type.
**Next steps:**
- List all unique `PTM.ModificationTitle` values and select PTM classes relevant to inflammasome regulation for first plots.
- Create `POI.csv` (UniProt accessions) for NLRP3/inflammasome pathway components.
- Build `scripts/01_load_and_triage.R` and `scripts/02_volcano_plots.R` and save outputs to `results/figures/`.

---

## LN-010 — Template for future entries

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

