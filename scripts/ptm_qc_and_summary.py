from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main() -> None:
    # --- Paths (assumes you run from repo root) ---
    repo_root = Path(".").resolve()
    infile = repo_root / "data_raw" / "PTM1.tsv"
    out_table_dir = repo_root / "tables"
    out_fig_dir = repo_root / "figures"

    out_table_dir.mkdir(parents=True, exist_ok=True)
    out_fig_dir.mkdir(parents=True, exist_ok=True)

    if not infile.exists():
        raise FileNotFoundError(
            f"Input file not found: {infile}\n"
            "Put PTM1.tsv in: data_raw/PTM1.tsv (relative to repo root)."
        )

    # --- Load ---
    df = pd.read_csv(infile, sep="\t", low_memory=False)

    # --- QC: shape ---
    print("=== QC SUMMARY ===")
    print(f"File: {infile}")
    print(f"Shape: {df.shape[0]:,} rows × {df.shape[1]:,} cols")
    print()

    # --- QC: missingness by column ---
    missing_pct = (df.isna().mean() * 100).sort_values(ascending=False)

    print("Top 20 columns by missingness (%):")
    print(missing_pct.head(20).to_string())
    print()

    missing_out = out_table_dir / "PTM1_missingness_by_column.csv"
    missing_pct.rename("missing_percent").to_csv(missing_out, index=True)
    print(f"Saved full missingness table -> {missing_out}")
    print()

    # --- QC: Qvalue distribution ---
    if "Qvalue" not in df.columns:
        raise KeyError("Expected column 'Qvalue' not found in PTM1.tsv")

    qv = pd.to_numeric(df["Qvalue"], errors="coerce")
    qv_valid = qv.dropna()

    print(f"Qvalue non-null: {qv_valid.shape[0]:,} / {df.shape[0]:,}")
    print("Qvalue quantiles:")
    print(qv_valid.quantile([0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]).to_string())
    print()

    plt.figure()
    plt.hist(qv_valid.clip(0, 1), bins=50)
    plt.xlabel("Qvalue (clipped to [0,1])")
    plt.ylabel("Count")
    plt.title("Qvalue distribution")
    qvalue_fig = out_fig_dir / "PTM1_qvalue_distribution.png"
    plt.tight_layout()
    plt.savefig(qvalue_fig, dpi=200)
    plt.close()
    print(f"Saved Qvalue histogram -> {qvalue_fig}")
    print()

    # Optional: -log10(Qvalue) view (useful if many tiny q-values)
    eps = 1e-300
    mlogq = -np.log10(qv_valid.clip(lower=eps))
    plt.figure()
    plt.hist(mlogq, bins=50)
    plt.xlabel("-log10(Qvalue)")
    plt.ylabel("Count")
    plt.title("-log10(Qvalue) distribution")
    mlogq_fig = out_fig_dir / "PTM1_minuslog10_qvalue_distribution.png"
    plt.tight_layout()
    plt.savefig(mlogq_fig, dpi=200)
    plt.close()
    print(f"Saved -log10(Qvalue) histogram -> {mlogq_fig}")
    print()

    # --- QC: '# of Ratios' distribution ---
    ratio_col = "# of Ratios"
    if ratio_col not in df.columns:
        raise KeyError(f"Expected column '{ratio_col}' not found in PTM1.tsv")

    nrat = pd.to_numeric(df[ratio_col], errors="coerce").dropna()

    print(f"'{ratio_col}' non-null: {nrat.shape[0]:,} / {df.shape[0]:,}")
    print(f"'{ratio_col}' quantiles:")
    print(nrat.quantile([0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]).to_string())
    print()

    plt.figure()
    bins = np.arange(0, int(nrat.max()) + 2) - 0.5
    plt.hist(nrat, bins=bins)
    plt.xlabel("# of Ratios")
    plt.ylabel("Count")
    plt.title("# of Ratios distribution")
    ratios_fig = out_fig_dir / "PTM1_num_ratios_distribution.png"
    plt.tight_layout()
    plt.savefig(ratios_fig, dpi=200)
    plt.close()
    print(f"Saved # of Ratios histogram -> {ratios_fig}")
    print()

    # --- PTM summary per Comparison × PTM type ---
    required = ["Comparison (group1/group2)", "PTM.ModificationTitle", "AVG Log2 Ratio", "Qvalue"]
    missing_required = [c for c in required if c not in df.columns]
    if missing_required:
        raise KeyError(f"Missing required columns for summary: {missing_required}")

    tmp = df.copy()
    tmp["AVG Log2 Ratio"] = pd.to_numeric(tmp["AVG Log2 Ratio"], errors="coerce")
    tmp["Qvalue"] = pd.to_numeric(tmp["Qvalue"], errors="coerce")

    alpha = 0.05
    tmp["_is_sig"] = tmp["Qvalue"].notna() & (tmp["Qvalue"] <= alpha)
    tmp["_is_up_sig"] = tmp["_is_sig"] & (tmp["AVG Log2 Ratio"] > 0)
    tmp["_is_down_sig"] = tmp["_is_sig"] & (tmp["AVG Log2 Ratio"] < 0)

    group_cols = ["Comparison (group1/group2)", "PTM.ModificationTitle"]

    summary = (
        tmp.groupby(group_cols, dropna=False)
        .agg(
            n=("Qvalue", "size"),
            n_qvalue_nonnull=("Qvalue", lambda s: s.notna().sum()),
            n_sig=("_is_sig", "sum"),
            n_up_sig=("_is_up_sig", "sum"),
            n_down_sig=("_is_down_sig", "sum"),
            median_abs_log2=("AVG Log2 Ratio", lambda s: np.nanmedian(np.abs(s.to_numpy()))),
            median_log2=("AVG Log2 Ratio", "median"),
            mean_log2=("AVG Log2 Ratio", "mean"),
            median_qvalue=("Qvalue", "median"),
        )
        .reset_index()
        .sort_values(["Comparison (group1/group2)", "PTM.ModificationTitle"])
        .reset_index(drop=True)
    )

    out_csv = out_table_dir / "PTM1_summary_by_comparison_ptmtype.csv"
    summary.to_csv(out_csv, index=False)
    print(f"Saved curated PTM summary table -> {out_csv}")
    print()
    print("Preview (first 20 rows):")
    print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
