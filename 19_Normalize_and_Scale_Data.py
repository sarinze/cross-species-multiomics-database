import pandas as pd
import numpy as np

INPUT_FILE = "Master_Multiomics_Matrix_With_KPPR1_TnSeq.csv"

OUTPUT_FILE = "Master_Multiomics_Matrix_Normalized.csv"
REPORT_FILE = "Normalization_Report.csv"
COLUMNS_FILE = "Normalized_Columns.txt"

print("\n========== STEP 19: NORMALIZE AND SCALE DATA ==========")

df = pd.read_csv(INPUT_FILE, low_memory=False)

print(f"\nLoaded file: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# Words that suggest a column SHOULD be normalized
include_keywords = [
    "log2fc",
    "log2_fold_change",
    "fold_change",
    "fitness",
    "effect_size",
    "selection_ratio",
    "mean_si",
    "rep1_si",
    "rep2_si",
    "rep3_si"
]

# Words that suggest a column should NOT be normalized
exclude_keywords = [
    "p_value",
    "pvalue",
    "adjusted_p",
    "padj",
    "q_value",
    "basemean",
    "count",
    "mean_input",
    "mean_output",
    "insertion",
    "length",
    "start",
    "end",
    "coordinate",
    "accession",
    "locus",
    "tag",
    "gene_name",
    "protein_name",
    "chromosome",
    "strand",
    "orientation",
    "annotation",
    "description",
    "product",
    "category",
    "type"
]

normalized_columns = []
report = []

df_norm = df.copy()

for col in df.columns:
    col_lower = col.lower()

    should_include = any(key in col_lower for key in include_keywords)
    should_exclude = any(key in col_lower for key in exclude_keywords)

    if should_include and not should_exclude:
        numeric_col = pd.to_numeric(df[col], errors="coerce")

        non_missing = numeric_col.notna().sum()

        if non_missing < 2:
            continue

        mean = numeric_col.mean()
        std = numeric_col.std()

        if std == 0 or pd.isna(std):
            continue

        z_scaled = (numeric_col - mean) / std

        df_norm[col] = z_scaled
        normalized_columns.append(col)

        report.append({
            "Column": col,
            "Non_Missing_Values": non_missing,
            "Original_Mean": mean,
            "Original_Std": std,
            "Normalized_Mean": z_scaled.mean(),
            "Normalized_Std": z_scaled.std(),
            "Missing_Values": numeric_col.isna().sum()
        })

df_norm.to_csv(OUTPUT_FILE, index=False)
pd.DataFrame(report).to_csv(REPORT_FILE, index=False)

with open(COLUMNS_FILE, "w") as f:
    for col in normalized_columns:
        f.write(col + "\n")

print("\n========== NORMALIZATION COMPLETE ==========")
print(f"Normalized columns: {len(normalized_columns)}")

print("\nColumns normalized:")
for col in normalized_columns:
    print(f"- {col}")

print(f"\nSaved normalized matrix: {OUTPUT_FILE}")
print(f"Saved normalization report: {REPORT_FILE}")
print(f"Saved normalized column list: {COLUMNS_FILE}")

print("\nImportant:")
print("Raw file was not changed.")
print("Missing values were preserved as NaN.")
print("P-values, BaseMean, IDs, annotations, coordinates, and protein lengths were not scaled.")