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

Change to Excel spreadsheet, values were off due to counting repeated evidence rows, thats why I was seeing multiples of 15. count of 15 = one site, being detected across different conditions, I think its actually being detected once in every condition, Not entirely sure what that means.
Decision: manipulate Data in R to create more accurate table, current table is misleading and uninformative.

Account for protein abundance changes in PTMs - See if there is actually a change in the modification or if it is just due to the abundance of the protein changing e.g IL1B protein abundance will increase in activated states -> shows increase in ptm but could just be increase of protein

change visual of my plots: need to decide this, when plotting modification/ comparison specific volcano plot, is it better to just plot those modifications, or plot all modifications but highlight them in a different colour with labels?

change from alphamap GUI to command line

exports from spectronaut, one including protein abundance change, and one with mass shift - check these against unimod to see if bioloigcal or artefact

genemania: upload a list of UniProt IDs to see their bioloigcal role and interactions with eachother.

