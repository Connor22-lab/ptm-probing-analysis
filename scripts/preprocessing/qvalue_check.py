import pandas as pd

df = pd.read_csv("data_raw/PTM1.tsv", sep="\t", low_memory=False)

print("Rows, Cols:", df.shape)

print("\nColumns containing 'q' or 'p':")
for c in df.columns:
    if "q" in c.lower() or "p" in c.lower():
        print(" -", c)

# -----------------------------
# Qvalue checks
# -----------------------------
q = pd.to_numeric(df["Qvalue"], errors="coerce")

print("\nQvalue checks:")
print("Non-null:", int(q.notna().sum()))
print("Min/Median/Max:", float(q.min()), float(q.median()), float(q.max()))
print("Unique Qvalues:", int(q.nunique()))

print("\nTop repeated Qvalues:")
print(q.value_counts(dropna=True).head(10).to_string())

# How many pass common thresholds (overall)
thresholds = [0.01, 0.05, 0.10, 0.20]
print("\nOverall counts by Qvalue threshold:")
for a in thresholds:
    n = int((q <= a).sum())
    pct = (n / len(q)) * 100 if len(q) else 0
    print(f"  Q <= {a:0.2f}: {n:,} ({pct:0.4f}%)")

# Per-comparison: where the signal is concentrated (if anywhere)
comp_col = "Comparison (group1/group2)"
if comp_col in df.columns:
    tmp = df[[comp_col]].copy()
    tmp["Qvalue"] = q

    for a in [0.05, 0.10]:
        per_comp = (
            tmp.groupby(comp_col)["Qvalue"]
               .apply(lambda s: int((s <= a).sum()))
               .sort_values(ascending=False)
        )
        print(f"\nTop comparisons by # of rows with Q <= {a:0.2f}:")
        print(per_comp.head(15).to_string())

# -----------------------------
# Pvalue checks (helps explain why Q is high)
# -----------------------------
if "Pvalue" in df.columns:
    p = pd.to_numeric(df["Pvalue"], errors="coerce")
    print("\nPvalue checks:")
    print("Non-null:", int(p.notna().sum()))
    print("Min/Median/Max:", float(p.min()), float(p.median()), float(p.max()))
    print("P <= 0.05:", int((p <= 0.05).sum()))
    print("P <= 0.01:", int((p <= 0.01).sum()))
