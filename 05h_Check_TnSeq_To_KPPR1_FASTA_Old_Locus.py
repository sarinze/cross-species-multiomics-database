import pandas as pd
import sqlite3

conn = sqlite3.connect("pseudomonas_new.db")

tnseq = pd.read_sql_query("""
SELECT Gene_ID, Gene_Name
FROM Klebsiella_TnSeq
""", conn)

# This file may contain old VK055 IDs from your earlier FASTA work
fasta_map = pd.read_csv("KPPR1_Protein_FASTA_With_Old_Locus.csv")

print("FASTA map columns:")
print(fasta_map.columns.tolist())

print("\nFirst 5 rows:")
print(fasta_map.head())

print("\nTnSeq examples:")
print(tnseq.head())

conn.close()