import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"
MATRIX_FILE = "Master_Matrix_With_Klebsiella_DEG.csv"

conn = sqlite3.connect(DB_NAME)

matrix = pd.read_csv(MATRIX_FILE)

tnseq = pd.read_sql_query("""
SELECT Gene_ID
FROM Klebsiella_TnSeq
WHERE Gene_ID IS NOT NULL
""", conn)

print("Examples from master linker Klebsiella_Locus_Tag:")
print(matrix["Klebsiella_Locus_Tag"].dropna().drop_duplicates().head(20).tolist())

print("\nExamples from Klebsiella_TnSeq Gene_ID:")
print(tnseq["Gene_ID"].dropna().drop_duplicates().head(20).tolist())

matrix_ids = set(matrix["Klebsiella_Locus_Tag"].dropna().astype(str))
tnseq_ids = set(tnseq["Gene_ID"].dropna().astype(str))

overlap = matrix_ids.intersection(tnseq_ids)

print("\nNumber of unique Klebsiella IDs in matrix:", len(matrix_ids))
print("Number of unique Gene_IDs in TnSeq:", len(tnseq_ids))
print("Number of overlapping IDs:", len(overlap))

if len(overlap) > 0:
    print("\nExample overlaps:")
    print(list(overlap)[:20])
else:
    print("\nNo overlap found.")

conn.close()