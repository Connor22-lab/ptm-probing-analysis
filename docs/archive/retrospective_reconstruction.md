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