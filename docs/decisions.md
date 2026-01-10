## Global decisions (living list)

1) **Primary analytical framing:** condition-based contrasts that reflect inflammasome biology (e.g., priming vs activation) rather than timepoints (unless timepoints are present in a future dataset).
2) **Initial significance thresholds (for exploratory volcano plots):**
   - log2FC cutoff: ±0.58 (≈1.5-fold)
   - Q-value cutoff: 0.05  
   These may be tightened later depending on PTM class, localisation confidence, and replication.
3) **POI labelling strategy:** proteins of interest stored in `POI.csv` as UniProt accessions; matching must handle semicolon-delimited `UniProtIds` values from Spectronaut.
4) **Reproducibility:** all plots saved with informative filenames including comparison + PTM type + thresholds.
5) Primary PTM class to focus next: Phospho (STY). Rationale: highest candidate density (rows_candidate and candidate_frac) among signalling-relevant PTMs under evidence filtering, with strong median effect size.

Secondary PTM class (follow-up): GlyGly (K). Rationale: larger median effect size but lower candidate fraction; may require more careful handling due to sparsity.