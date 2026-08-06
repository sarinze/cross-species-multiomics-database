import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"
INPUT_FILE = "Master_Functional_Genomics_Matrix.csv"
OUTPUT_FILE = "Master_Matrix_With_Klebsiella_DEG.csv"

conn = sqlite3.connect(DB_NAME)

matrix = pd.read_csv(INPUT_FILE)

kleb_deg = pd.read_sql_query("""
SELECT
    Gene_ID AS Klebsiella_Locus_Tag,

    K56_vs_Colistin_Log2FC,
    K56_vs_Colistin_P_Value,
    K56_vs_Colistin_Q_Value,

    K56_vs_Combination_Log2FC,
    K56_vs_Combination_P_Value,
    K56_vs_Combination_Q_Value,

    Colistin_vs_Combination_Log2FC,
    Colistin_vs_Combination_P_Value,
    Colistin_vs_Combination_Q_Value
FROM KPHS_K56_COLISTIN_COMBINATION_DEG
""", conn)

kleb_deg = kleb_deg.drop_duplicates(
    subset=["Klebsiella_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    kleb_deg,
    on="Klebsiella_Locus_Tag",
    how="left"
)

print("After adding Klebsiella DEG:", matrix.shape)

matrix.to_csv(OUTPUT_FILE, index=False)

print("\nKlebsiella DEG added successfully")
print("Output file:", OUTPUT_FILE)
print("Rows:", len(matrix))
print("Columns:", len(matrix.columns))

print("\nMissing DEG values:")
print(matrix[
    [
        "K56_vs_Colistin_Log2FC",
        "K56_vs_Colistin_Q_Value",
        "K56_vs_Combination_Log2FC",
        "K56_vs_Combination_Q_Value"
    ]
].isna().sum())

conn.close()