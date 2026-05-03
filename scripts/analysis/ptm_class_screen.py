import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_TSV = PROJECT_ROOT / "data_raw" / "PTM1.tsv"
OUTDIR = PROJECT_ROOT / "results" / "summary"      

# === FILTER SETTINGS (tweak later if needed) ===
MIN_RATIOS = 2
MIN_UNIQUE_PEPTIDES = 2
MIN_ABS_LOG2 = 0.3
MAX_PVALUE = 0.05

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT_TSV, sep="\t", low_memory=False)

    # Coerce numeric columns
    num_cols = [
        "AVG Log2 Ratio", "Absolute AVG Log2 Ratio",
        "Pvalue", "Qvalue", "# of Ratios", "# Unique Total Peptides"
    ]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # --- Targeted shortlist: UniProt Based ---
    ptm_keep = ["Phospho (STY)", "GlyGly (K)", "Acetyl (K)", "HexNAc (ST)"]

    uniprot_keep = [
        "Q96P20",  # NLRP3
        "Q9ULZ3",  # PYCARD / ASC
        "P29466",  # CASP1
        "P01584",  # IL1B
        "Q14116",  # IL18
        "P57764",  # GSDMD
        "Q8TDX7",  # NEK7
        "O14862",  # AIM2
        "Q9NPP4",  # NLRC4
        "Q06187",  # BTK
    ]

    f = df[df["PTM.ModificationTitle"].isin(ptm_keep)].copy()
    pattern = "|".join(uniprot_keep)
    f = f[f["UniProtIds"].astype(str).str.contains(pattern, na=False)].copy()

    # Optional relaxed evidence filter so the shortlist isn't pure noise
    f = f[(f["# of Ratios"] >= MIN_RATIOS) & (f["Absolute AVG Log2 Ratio"] >= MIN_ABS_LOG2)].copy()

    f = f.sort_values(
        ["Pvalue", "Absolute AVG Log2 Ratio", "# of Ratios", "# Unique Total Peptides"],
        ascending=[True, False, False, False]
    )

    f.to_csv(OUTDIR / "filtered__targets_UniProt.csv", index=False)
    print("Wrote:", OUTDIR / "filtered__targets_UniProt.csv", "rows:", len(f))
    # --- Targeted shortlist (GENE-based): inflammasome genes + PTM classes of interest ---
    ptm_keep = ["Phospho (STY)", "GlyGly (K)", "Acetyl (K)", "HexNAc (ST)"]
    genes_keep = ["NLRP3", "PYCARD", "CASP1", "GSDMD", "IL1B", "IL18", "NEK7", "AIM2", "NLRC4", "BTK"]

    g = df[df["PTM.ModificationTitle"].isin(ptm_keep)].copy()
    g = g[g["Genes"].astype(str).str.contains("|".join(genes_keep), na=False)].copy()

    # light evidence filter
    g = g[(g["# of Ratios"] >= MIN_RATIOS) & (g["Absolute AVG Log2 Ratio"] >= MIN_ABS_LOG2)].copy()

    g = g.sort_values(
        ["Pvalue", "Absolute AVG Log2 Ratio", "# of Ratios", "# Unique Total Peptides"],
        ascending=[True, False, False, False]
    )

    g.to_csv(OUTDIR / "filtered_targets_genes.csv", index=False)
    print("Wrote:", OUTDIR / "filtered_targets_genes.csv", "rows:", len(g))

    # Evidence filter (power filter)
    e = df[
        (df["# of Ratios"] >= MIN_RATIOS) &
        (df["# Unique Total Peptides"] >= MIN_UNIQUE_PEPTIDES) &
        (df["Absolute AVG Log2 Ratio"] >= MIN_ABS_LOG2)
    ].copy()

    # Summary by PTM class
    summary = (
        e.groupby("PTM.ModificationTitle")
         .agg(
             rows_total=("ProteinGroups", "size"),
             rows_candidate=("Pvalue", lambda s: (s <= MAX_PVALUE).sum()),
             candidate_frac=("Pvalue", lambda s: (s <= MAX_PVALUE).mean()),
             median_abs_log2=("Absolute AVG Log2 Ratio", "median"),
             median_p=("Pvalue", "median"),
             median_ratios=("# of Ratios", "median"),
             median_unique_peps=("# Unique Total Peptides", "median"),
         )
         .sort_values(["rows_candidate", "candidate_frac", "median_abs_log2"], ascending=False)
    )

    # Summary by comparison + PTM class
    summary_by_comp = (
        e.groupby(["Comparison (group1/group2)", "PTM.ModificationTitle"])
         .agg(
             rows_total=("ProteinGroups", "size"),
             rows_candidate=("Pvalue", lambda s: (s <= MAX_PVALUE).sum()),
             candidate_frac=("Pvalue", lambda s: (s <= MAX_PVALUE).mean()),
             median_abs_log2=("Absolute AVG Log2 Ratio", "median"),
             median_p=("Pvalue", "median"),
             median_ratios=("# of Ratios", "median"),
         )
         .sort_values(["rows_candidate", "candidate_frac", "median_abs_log2"], ascending=False)
    )

    # Save outputs
    summary.to_csv(OUTDIR / "ptm_class_summary.csv")
    summary_by_comp.to_csv(OUTDIR / "ptm_class_summary_by_comparison.csv")

    print("Wrote:")
    print(" -", OUTDIR / "ptm_class_summary.csv")
    print(" -", OUTDIR / "ptm_class_summary_by_comparison.csv")
    print("\nTop PTM classes:")
    print(summary.head(10))

if __name__ == "__main__":
    main()
