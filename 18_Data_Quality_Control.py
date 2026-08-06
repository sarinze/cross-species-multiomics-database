import pandas as pd
from pathlib import Path

# =========================
# SETTINGS
# =========================

MASTER_FILE = "Master_Multiomics_Matrix_With_KPPR1_TnSeq.csv"

OUTPUT_FOLDER = Path("18_QC_Results")
OUTPUT_FOLDER.mkdir(exist_ok=True)

# =========================
# LOAD MASTER MATRIX
# =========================

print("\n========== DATA QUALITY CONTROL ==========")

df = pd.read_csv(MASTER_FILE)

print(f"\nLoaded file: {MASTER_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# =========================
# 1. BASIC STRUCTURE CHECK
# =========================

basic_summary = pd.DataFrame({
    "Column": df.columns,
    "Data_Type": df.dtypes.astype(str).values,
    "Non_Null_Count": df.notna().sum().values,
    "Null_Count": df.isna().sum().values,
    "Percent_Missing": (df.isna().mean() * 100).round(2).values,
    "Unique_Values": df.nunique(dropna=True).values
})

basic_summary.to_csv(OUTPUT_FOLDER / "QC_Column_Missingness_Summary.csv", index=False)

print("\nSaved: QC_Column_Missingness_Summary.csv")

# =========================
# 2. DUPLICATE FULL ROW CHECK
# =========================

duplicate_rows = df[df.duplicated()]

duplicate_rows.to_csv(OUTPUT_FOLDER / "QC_Duplicate_Full_Rows.csv", index=False)

print(f"\nDuplicate full rows: {duplicate_rows.shape[0]}")
print("Saved: QC_Duplicate_Full_Rows.csv")

# =========================
# 3. DUPLICATE GENE ID CHECKS
# =========================

possible_id_cols = [
    "PA14_Original_Locus_Tag",
    "PA14_New_Locus_Tag",
    "PAO1_Locus_Tag",
    "MGH78578_Locus_Tag",
    "NCTC5055_Locus_Tag",
    "HS11286_Locus_Tag",
    "KPPR1_Locus_Tag",
    "KPPR1_Gene_ID",
    "KPPR1_Gene"
]

existing_id_cols = [col for col in possible_id_cols if col in df.columns]

dup_gene_summary = []

for col in existing_id_cols:
    duplicated = df[df[col].notna() & df[col].duplicated(keep=False)]
    
    duplicated.to_csv(
        OUTPUT_FOLDER / f"QC_Duplicates_{col}.csv",
        index=False
    )

    dup_gene_summary.append({
        "Column": col,
        "Duplicated_Rows": duplicated.shape[0],
        "Duplicated_Unique_IDs": duplicated[col].nunique()
    })

dup_gene_summary_df = pd.DataFrame(dup_gene_summary)
dup_gene_summary_df.to_csv(
    OUTPUT_FOLDER / "QC_Duplicate_Gene_ID_Summary.csv",
    index=False
)

print("\nSaved duplicate gene ID checks.")

# =========================
# 4. NUMERIC COLUMN SUMMARY
# =========================

numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()

if numeric_cols:
    numeric_summary = df[numeric_cols].describe().T.reset_index()
    numeric_summary = numeric_summary.rename(columns={"index": "Column"})

    numeric_summary.to_csv(
        OUTPUT_FOLDER / "QC_Numeric_Column_Summary.csv",
        index=False
    )

    print("\nSaved: QC_Numeric_Column_Summary.csv")
else:
    print("\nNo numeric columns found.")

# =========================
# 5. OMICS VALUE CHECKS
# =========================

omics_keywords = [
    "Log2FC",
    "log2fc",
    "Fold_Change",
    "fold_change",
    "P_Value",
    "p_value",
    "Adjusted_P_Value",
    "adjusted",
    "padj",
    "BaseMean",
    "Fitness",
    "TnSeq"
]

omics_cols = [
    col for col in df.columns
    if any(keyword.lower() in col.lower() for keyword in omics_keywords)
]

omics_summary = []

for col in omics_cols:
    numeric = pd.to_numeric(df[col], errors="coerce")

    omics_summary.append({
        "Column": col,
        "Non_Null_Count": numeric.notna().sum(),
        "Null_Count": numeric.isna().sum(),
        "Percent_Missing": round(numeric.isna().mean() * 100, 2),
        "Minimum": numeric.min(),
        "Maximum": numeric.max(),
        "Mean": numeric.mean(),
        "Median": numeric.median()
    })

omics_summary_df = pd.DataFrame(omics_summary)
omics_summary_df.to_csv(
    OUTPUT_FOLDER / "QC_Omics_Value_Summary.csv",
    index=False
)

print("\nSaved: QC_Omics_Value_Summary.csv")

# =========================
# 6. NEAR EMPTY COLUMN CHECK
# =========================

near_empty_cols = basic_summary[basic_summary["Percent_Missing"] >= 95]

near_empty_cols.to_csv(
    OUTPUT_FOLDER / "QC_Near_Empty_Columns_95percent.csv",
    index=False
)

print(f"\nColumns with >=95% missing values: {near_empty_cols.shape[0]}")
print("Saved: QC_Near_Empty_Columns_95percent.csv")

# =========================
# 7. OMICS COVERAGE PER ROW
# =========================

if omics_cols:
    df["QC_Omics_Non_Null_Count"] = df[omics_cols].notna().sum(axis=1)

    coverage_cols = existing_id_cols + ["QC_Omics_Non_Null_Count"]

    row_coverage = df[coverage_cols]

    row_coverage.to_csv(
        OUTPUT_FOLDER / "QC_Row_Omics_Coverage.csv",
        index=False
    )

    print("\nSaved: QC_Row_Omics_Coverage.csv")
else:
    print("\nNo omics columns detected for row coverage.")

# =========================
# 8. ROWS WITH ZERO OMICS DATA
# =========================

if omics_cols:
    zero_omics_rows = df[df["QC_Omics_Non_Null_Count"] == 0]

    zero_omics_rows.to_csv(
        OUTPUT_FOLDER / "QC_Rows_With_Zero_Omics_Data.csv",
        index=False
    )

    print(f"\nRows with zero omics data: {zero_omics_rows.shape[0]}")
    print("Saved: QC_Rows_With_Zero_Omics_Data.csv")

# =========================
# 9. EXPORT QC MASTER COPY
# =========================

qc_output_csv = OUTPUT_FOLDER / "Master_Multiomics_Matrix_WITH_QC_COLUMNS.csv"

df.to_csv(qc_output_csv, index=False)

print(f"\nSaved QC master copy: {qc_output_csv}")

print("\n========== QC COMPLETE ==========")
print(f"All QC files saved in: {OUTPUT_FOLDER}")