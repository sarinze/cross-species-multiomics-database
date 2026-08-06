import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"

conn = sqlite3.connect(DB_NAME)

tnseq = pd.read_sql_query("""
SELECT
    Gene_ID,
    Gene_Name,
    Gene_Length,
    Total_Insertions,
    Unique_Insertions,
    Mean_Input,
    Mean_Output,
    Fold_Change,
    Log2_Fold_Change_Output_Input,
    P_Value,
    Adjusted_P_Value
FROM Klebsiella_TnSeq
""", conn)

print("Total rows:", len(tnseq))
print("Unique Gene_IDs:", tnseq["Gene_ID"].nunique())

gene_counts = tnseq["Gene_ID"].value_counts()

print("\nRows per gene summary:")
print(gene_counts.describe())

print("\nTop 20 genes with most rows:")
print(gene_counts.head(20))

print("\nFirst 10 rows:")
print(tnseq.head(10))

conn.close()