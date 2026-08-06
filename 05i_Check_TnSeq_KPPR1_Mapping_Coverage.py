import pandas as pd
import sqlite3

conn = sqlite3.connect("pseudomonas_new.db")

tnseq = pd.read_sql_query("""
SELECT *
FROM Klebsiella_TnSeq
""", conn)

fasta_map = pd.read_csv("KPPR1_Protein_FASTA_With_Old_Locus.csv")

# Keep only rows with old locus tags
fasta_map = fasta_map.dropna(subset=["Old_Locus_Tag"])

# Avoid duplicate old locus tags
fasta_map = fasta_map.drop_duplicates(
    subset=["Old_Locus_Tag"],
    keep="first"
)

merged = tnseq.merge(
    fasta_map[["Old_Locus_Tag", "Protein_Accession"]],
    left_on="Gene_ID",
    right_on="Old_Locus_Tag",
    how="left"
)

print("TnSeq rows:", len(tnseq))
print("Mapped to KPPR1 protein accession:", merged["Protein_Accession"].notna().sum())
print("Unmapped:", merged["Protein_Accession"].isna().sum())

print("\nFirst 10 mapped rows:")
print(merged[merged["Protein_Accession"].notna()].head(10)[
    ["Gene_ID", "Gene_Name", "Protein_Accession"]
])

print("\nFirst 10 unmapped rows:")
print(merged[merged["Protein_Accession"].isna()].head(10)[
    ["Gene_ID", "Gene_Name"]
])

conn.close()