import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql_query("""
SELECT
    Locus_Tag,
    Gene_Type,
    Protein_Name,
    Protein_Accession,
    Protein_Length,
    FASTA_Protein_Name,
    FASTA_Protein_Length
FROM kppr1_annotation_with_refseq_fasta
WHERE Protein_Sequence IS NULL
""", conn)

conn.close()

print("Missing rows:", len(df))
print("\nGene type counts:")
print(df["Gene_Type"].value_counts(dropna=False))

print("\nFirst 50 missing rows:")
print(df.head(50))

df.to_csv("KPPR1_Missing_FASTA_Sequences.csv", index=False)

print("\nSaved CSV: KPPR1_Missing_FASTA_Sequences.csv")