import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

df = pd.read_sql("SELECT * FROM PAO1_PA14_RBH_LINKER", conn)

clean = df.rename(columns={
    "PAO1_PROTEIN": "PAO1_SEQUENCE",
    "PA14_PROTEIN": "PA14_SEQUENCE",

    "Locus_Tag": "PAO1_Locus_Tag",
    "Original_Locus_Tag": "PA14_Original_Locus_Tag",
    "New_Locus_Tag": "PA14_New_Locus_Tag",

    "Genome_Accession_PAO1": "PAO1_Genome_Accession",
    "Start_PAO1": "PAO1_Start",
    "End_PAO1": "PAO1_End",
    "Orientation_PAO1": "PAO1_Orientation",
    "Chromosome_PAO1": "PAO1_Chromosome",
    "Gene_Name_PAO1": "PAO1_Gene_Name",
    "Symbol_PAO1": "PAO1_Symbol",
    "Gene_Type_PAO1": "PAO1_Gene_Type",
    "Protein_Name_PAO1": "PAO1_Protein_Name",
    "Protein_Accession_PAO1": "PAO1_Protein_Accession",
    "Protein_Length_PAO1": "PAO1_Protein_Length",

    "Genome_Accession_PA14": "PA14_Genome_Accession",
    "Start_PA14": "PA14_Start",
    "End_PA14": "PA14_End",
    "Orientation_PA14": "PA14_Orientation",
    "Chromosome_PA14": "PA14_Chromosome",
    "Gene_Name_PA14": "PA14_Gene_Name",
    "Symbol_PA14": "PA14_Symbol",
    "Gene_Type_PA14": "PA14_Gene_Type",
    "Protein_Name_PA14": "PA14_Protein_Name",
    "Protein_Accession_PA14": "PA14_Protein_Accession",
    "Protein_Length_PA14": "PA14_Protein_Length",
})

clean = clean[[
    "PA14_SEQUENCE",
    "PA14_Original_Locus_Tag",
    "PA14_New_Locus_Tag",
    "PA14_Gene_Name",
    "PA14_Symbol",
    "PA14_Gene_Type",
    "PA14_Protein_Name",
    "PA14_Protein_Accession",
    "PA14_Protein_Length",
    "PA14_Genome_Accession",
    "PA14_Start",
    "PA14_End",
    "PA14_Orientation",
    "PA14_Chromosome",

    "PAO1_SEQUENCE",
    "PAO1_Locus_Tag",
    "PAO1_Gene_Name",
    "PAO1_Symbol",
    "PAO1_Gene_Type",
    "PAO1_Protein_Name",
    "PAO1_Protein_Accession",
    "PAO1_Protein_Length",
    "PAO1_Genome_Accession",
    "PAO1_Start",
    "PAO1_End",
    "PAO1_Orientation",
    "PAO1_Chromosome",

    "PERCENT_IDENTITY",
    "BITSCORE"
]]

clean.to_csv("PA14_PAO1_CLEAN_LINKER.csv", index=False)

clean.to_sql(
    "PA14_PAO1_CLEAN_LINKER",
    conn,
    if_exists="replace",
    index=False
)

print("Saved CSV: PA14_PAO1_CLEAN_LINKER.csv")
print("Imported into SQLite as: PA14_PAO1_CLEAN_LINKER")
print("Rows:", len(clean))
print("Columns:")
print(clean.columns.tolist())

conn.close()