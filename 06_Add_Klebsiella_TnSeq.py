import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"
INPUT_FILE = "Master_Matrix_With_Klebsiella_DEG.csv"
OUTPUT_FILE = "Master_Multiomics_Matrix.csv"

conn = sqlite3.connect(DB_NAME)

matrix = pd.read_csv(INPUT_FILE)

kleb_tnseq = pd.read_sql_query("""
SELECT
    Gene_ID AS Klebsiella_Locus_Tag,
    Total_Insertions AS Klebsiella_TnSeq_Total_Insertions,
    Unique_Insertions AS Klebsiella_TnSeq_Unique_Insertions,
    Mean_Input AS Klebsiella_TnSeq_Mean_Input,
    Mean_Output AS Klebsiella_TnSeq_Mean_Output,
    Fold_Change AS Klebsiella_TnSeq_Fold_Change,
    Log2_Fold_Change_Output_Input AS Klebsiella_TnSeq_Log2FC_Output_Input,
    P_Value AS Klebsiella_TnSeq_P_Value,
    Adjusted_P_Value AS Klebsiella_TnSeq_Adjusted_P_Value
FROM Klebsiella_TnSeq
""", conn)

kleb_tnseq = kleb_tnseq.drop_duplicates(
    subset=["Klebsiella_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    kleb_tnseq,
    on="Klebsiella_Locus_Tag",
    how="left"
)

print("After adding Klebsiella Tn-seq:", matrix.shape)

matrix.to_csv(OUTPUT_FILE, index=False)

print("\nMaster multi-omics matrix created successfully")
print("Output file:", OUTPUT_FILE)
print("Rows:", len(matrix))
print("Columns:", len(matrix.columns))

print("\nMissing Klebsiella Tn-seq values:")
print(matrix[
    [
        "Klebsiella_TnSeq_Log2FC_Output_Input",
        "Klebsiella_TnSeq_Adjusted_P_Value"
    ]
].isna().sum())

conn.close()