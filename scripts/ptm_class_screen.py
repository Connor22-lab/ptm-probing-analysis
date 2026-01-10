import pandas as pd
from pathlib import Path

# === EDIT THESE TWO LINES ===
INPUT_TSV = r"C:\Users\conno\Desktop\Final Year School\Final Year Project\data_raw\PTM1.tsv"   
OUTDIR = Path("results_small")      

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
