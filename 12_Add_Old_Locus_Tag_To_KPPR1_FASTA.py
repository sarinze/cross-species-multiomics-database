import sqlite3
import pandas as pd
import re

DB_FILE = "pseudomonas_new.db"

conn = sqlite3.connect(DB_FILE)

df = pd.read_sql_query("""
SELECT *
FROM kppr1_protein_fasta
""", conn)

def extract_old_locus(header):
    match = re.search(r"VK055_\d+", str(header))
    return match.group(0) if match else None

df["Old_Locus_Tag"] = df["Full_FASTA_Header"].apply(extract_old_locus)

print("Rows:", len(df))
print("Rows with old locus tag:", df["Old_Locus_Tag"].notna().sum())
print("Rows without old locus tag:", df["Old_Locus_Tag"].isna().sum())

print("\nFirst 30 rows:")
print(df[["Protein_Accession", "Old_Locus_Tag", "Protein_Name"]].head(30))

df.to_csv("KPPR1_Protein_FASTA_With_Old_Locus.csv", index=False)

df.to_sql(
    "kppr1_protein_fasta_with_old_locus",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nSaved CSV: KPPR1_Protein_FASTA_With_Old_Locus.csv")
print("Imported SQLite table: kppr1_protein_fasta_with_old_locus")