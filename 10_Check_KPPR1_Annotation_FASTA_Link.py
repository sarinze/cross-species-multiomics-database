import sqlite3
import pandas as pd

DB_FILE = "pseudomonas_new.db"

conn = sqlite3.connect(DB_FILE)

query = """
SELECT
    a.*,
    f.Protein_Name AS FASTA_Protein_Name,
    f.Protein_Length AS FASTA_Protein_Length,
    f.Protein_Sequence,
    f.Full_FASTA_Header
FROM KPPR1_Genome_Annotation a
LEFT JOIN kppr1_refseq_protein_fasta f
    ON a.Protein_Accession = f.Protein_Accession
"""

df = pd.read_sql_query(query, conn)

print("Rows after join:", len(df))
print("Matched FASTA sequences:", df["Protein_Sequence"].notna().sum())
print("Missing FASTA sequences:", df["Protein_Sequence"].isna().sum())

print("\nFirst 20 joined rows:")
print(df[[
    "Locus_Tag",
    "Gene_Name",
    "Protein_Accession",
    "Protein_Name",
    "FASTA_Protein_Name",
    "Protein_Length",
    "FASTA_Protein_Length"
]].head(20))

df.to_csv("KPPR1_Annotation_With_RefSeq_FASTA.csv", index=False)

df.to_sql(
    "kppr1_annotation_with_refseq_fasta",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nSaved CSV: KPPR1_Annotation_With_RefSeq_FASTA.csv")
print("Imported SQLite table: kppr1_annotation_with_refseq_fasta")