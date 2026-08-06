import pandas as pd
import numpy as np

print("\n========== STEP 23b: FILTER INFORMATIVE CLUSTERING INPUT ==========\n")

INPUT_FILE = "Clustering_Input_Matrix.csv"
OUTPUT_FILE = "Clustering_Input_Matrix_Filtered.csv"
FILTER_REPORT_FILE = "Clustering_Input_Filtering_Report.csv"

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows before filtering: {df.shape[0]}")
print(f"Columns before filtering: {df.shape[1]}")

# -----------------------------
# Identify ID / metadata columns
# -----------------------------
metadata_keywords = [
    "gene",
    "locus",
    "tag",
    "name",
    "description",
    "product",
    "protein",
    "accession",
    "cluster",
]

metadata_cols = [
    col for col in df.columns
    if any(keyword in col.lower() for keyword in metadata_keywords)
]

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

omics_cols = [
    col for col in numeric_cols
    if col not in metadata_cols
    and col != "Measured_Omics_Count"
]

print(f"\nDetected metadata columns: {len(metadata_cols)}")
print(f"Detected numeric omics columns: {len(omics_cols)}")

# -----------------------------
# Treat 0 as uninformative if matrix was zero-filled
# -----------------------------
omics_df = df[omics_cols].copy()

real_value_mask = omics_df.notna() & (omics_df != 0)

df["Measured_Omics_Count"] = real_value_mask.sum(axis=1)

# -----------------------------
# Row filtering
# -----------------------------
MIN_MEASURED_OMICS_PER_GENE = 6

filtered_df = df[df["Measured_Omics_Count"] >= MIN_MEASURED_OMICS_PER_GENE].copy()

print(f"\nMinimum measured omics values per gene: {MIN_MEASURED_OMICS_PER_GENE}")
print(f"Rows after row filtering: {filtered_df.shape[0]}")
print(f"Rows removed: {df.shape[0] - filtered_df.shape[0]}")

# -----------------------------
# Column filtering
# -----------------------------
MIN_NONZERO_FRACTION_PER_COLUMN = 0.05  # keep columns with at least 1% informative values

column_report = []

cols_to_keep = []

for col in omics_cols:
    nonzero_count = (filtered_df[col].notna() & (filtered_df[col] != 0)).sum()
    nonzero_fraction = nonzero_count / len(filtered_df) if len(filtered_df) > 0 else 0

    keep = nonzero_fraction >= MIN_NONZERO_FRACTION_PER_COLUMN

    column_report.append({
        "Column": col,
        "Nonzero_Count_After_Row_Filtering": nonzero_count,
        "Nonzero_Fraction_After_Row_Filtering": nonzero_fraction,
        "Kept": keep
    })

    if keep:
        cols_to_keep.append(col)

column_report_df = pd.DataFrame(column_report)

# Keep metadata + filtered omics + measured count
final_cols = metadata_cols + cols_to_keep + ["Measured_Omics_Count"]

# Remove duplicate column names if any
final_cols = list(dict.fromkeys([col for col in final_cols if col in filtered_df.columns]))

final_df = filtered_df[final_cols].copy()

print(f"\nOmics columns kept: {len(cols_to_keep)}")
print(f"Omics columns removed: {len(omics_cols) - len(cols_to_keep)}")
print(f"Final rows: {final_df.shape[0]}")
print(f"Final columns: {final_df.shape[1]}")

# -----------------------------
# Save outputs
# -----------------------------
final_df.to_csv(OUTPUT_FILE, index=False)
column_report_df.to_csv(FILTER_REPORT_FILE, index=False)

print(f"\nSaved filtered clustering input: {OUTPUT_FILE}")
print(f"Saved filtering report: {FILTER_REPORT_FILE}")

# -----------------------------
# Preview
# -----------------------------
print("\nMeasured omics count summary:")
print(final_df["Measured_Omics_Count"].describe())

print("\nTop retained omics columns:")
print(cols_to_keep[:20])

print("\n========== STEP 23b COMPLETE ==========\n")