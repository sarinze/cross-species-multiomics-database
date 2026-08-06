import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"
OUTPUT_FILE = "Master_Ortholog_Backbone.csv"

conn = sqlite3.connect(DB_NAME)

backbone = pd.read_sql_query("""
SELECT
    Pseudomonas_Strain,
    Klebsiella_Strain,
    Pseudomonas_Locus_Tag,
    Pseudomonas_Original_Locus_Tag,
    Pseudomonas_Gene_Name,
    Pseudomonas_Symbol,
    Pseudomonas_Protein_Accession,
    Klebsiella_Locus_Tag,
    Klebsiella_Gene_Name,
    Klebsiella_Symbol,
    Klebsiella_Protein_Accession,
    Percent_Identity,
    Alignment_Length,
    Evalue,
    Bitscore,
    Source_File
FROM KLEBSIELLA_PSEUDOMONAS_MASTER_LINKER
""", conn)

# Remove exact duplicate rows
backbone = backbone.drop_duplicates()

# Create unique ID for each Pseudomonas-Klebsiella ortholog pair
backbone["Ortholog_Pair_ID"] = (
    backbone["Pseudomonas_Strain"].astype(str) + "_" +
    backbone["Pseudomonas_Locus_Tag"].astype(str) + "__" +
    backbone["Klebsiella_Strain"].astype(str) + "_" +
    backbone["Klebsiella_Locus_Tag"].astype(str)
)

# Move ID column to the front
cols = ["Ortholog_Pair_ID"] + [c for c in backbone.columns if c != "Ortholog_Pair_ID"]
backbone = backbone[cols]

# Save output
backbone.to_csv(OUTPUT_FILE, index=False)

print("Master ortholog backbone created successfully")
print("Output file:", OUTPUT_FILE)
print("Rows:", len(backbone))
print("Columns:", len(backbone.columns))
print(backbone.head())

conn.close()