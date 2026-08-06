import pandas as pd
import sqlite3

csv_file = "KPPR1_Genome_Annotation.csv"
db_file = "pseudomonas_new.db"
table_name = "KPPR1_Genome_Annotation"

df = pd.read_csv(csv_file)

# Remove fully empty rows just in case
df = df.dropna(how="all")

# Standardize expected columns
df.columns = [
    "Locus_Tag",
    "Genome_Accession",
    "Start",
    "End",
    "Orientation",
    "Chromosome",
    "Gene_Name",
    "Symbol",
    "Gene_Type",
    "Protein_Name",
    "Protein_Accession",
    "Protein_Length"
]

print("========== QUALITY CONTROL ==========")
print("Rows:", len(df))
print("Columns:", len(df.columns))
print("Missing Locus_Tag:", df["Locus_Tag"].isna().sum())
print("Duplicate Locus_Tag:", df["Locus_Tag"].duplicated().sum())
print("Missing Protein_Accession:", df["Protein_Accession"].isna().sum())

conn = sqlite3.connect(db_file)

df.to_sql(
    table_name,
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\n========== IMPORT COMPLETE ==========")
print("Table:", table_name)
print("Rows imported:", len(df))
print("Columns:", list(df.columns))