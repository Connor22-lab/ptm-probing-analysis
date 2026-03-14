Workflow for intermediate anaylsis.

So have categorised PTM modifications, Unique proteins and residues, now we want to look at them further.
Main things we are looking for: Modifications with high levels of site localsiation confidence.
Modifications that change with high abundance across conditions.
Modifications that still change with high abundance when taking protein abundance into account.

So we have 3 files: PTM1.tsv: Our ptm candidates, site location identifier is group. these are a comparison conditions
ptm_open_search_site_localisation.tsv: this is a run wise with replicates, gives us site location and localisation score, will not be able to merge these two files, so will look at ptm_open_search_site_localisation first
ptm_open_search_PG_candidates: this is our protein candidates list, we will be merging this with PTM1.tsv to account for differential abundance.

So first all 3 files were loaded into R and their columns were printed, to see if I could merge PTM1.tsv with the site localsiation, because of the descripency between comparison between 2 groups no.


what I done next was filter all contaminants out of the dataset (Cont) in the group key
Then filtered out all artifacts.
 Retained 11 modification types: Phospho (STY), GlyGly (K), Acetyl (K), HexNAc (ST), Trimethyl (K), Crotonyl, Myristoyl, Hex (K), Sulfo (STY), Nitrosyl, and Oxidation (W).
Site localisation cut off was already applied in spectronaut.
Next step is to flip numerator and denominator as they do not reflect the actual biology + flip AVG log2 ratio in flipped comparisons so its consistent.
Then subset the dataset to 6 key biological comparison, unstim -> primed -> activated
Each PTM class quantified protein coverage (unique proteins) site coverage (unique sites) and mean absolute log2FC across these compairsons
Heat maps generated of protein coverage and effect size per PTM class per comparison
Volcano plots for each key compairson with inflammasome proteins of interest
POI focused analysis identifying which canonical and regulatory inflammasome proteins carried detectable modifications

Generated a fomral selection table ranking all 11 biological PTM classes by mean absolute log2FC across key comparisons, protein coverage, presence in all 6 comparisons and detection on inflammasome protein of interest.
Primary closed search candidates: Phospho, Acetyl, GlyGly, Trimethyl. selection on the basis of hjgih effect sizes, wide protein coverage, detection on canonical inflammasome proteins and then aswell strong prior biological evidence in the inflammasome literature.

The next step is to merge the ptm dataset with the protein candidates to calculate a corrected ptm log2fc, taking into account protein abundance.

