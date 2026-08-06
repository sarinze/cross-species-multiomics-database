import pandas as pd
import sqlite3

rbh = pd.read_csv("PAO1_NCTC5055_RBH_Orthologs.csv")
pao1 = pd.read_csv("PAO1_Genome_Annotation.csv")
nctc = pd.read_csv("NCTC5055_Genome_Annotation.csv")

pao1_keep = pao1[
    [
        "Locus_Tag",
        "Gene_Name",
        "Symbol",
        "Gene_Type",
        "Protein_Name",
        "Protein_Accession",
        "Protein_Length",
    ]
].copy()

pao1_keep = pao1_keep.rename(
    columns={
        "Locus_Tag": "PAO1_Locus_Tag",
        "Gene_Name": "PAO1_Gene_Name",
        "Symbol": "PAO1_Symbol",
        "Gene_Type": "PAO1_Gene_Type",
        "Protein_Name": "PAO1_Protein_Name",
        "Protein_Accession": "PAO1_Protein_Accession",
        "Protein_Length": "PAO1_Protein_Length",
    }
)

nctc_keep = nctc[
    [
        "Locus_Tag",
        "Gene_Name",
        "Symbol",
        "Gene_Type",
        "Protein_Name",
        "Protein_Accession",
        "Protein_Length",
    ]
].copy()

nctc_keep = nctc_keep.rename(
    columns={
        "Locus_Tag": "NCTC5055_Locus_Tag",
        "Gene_Name": "NCTC5055_Gene_Name",
        "Symbol": "NCTC5055_Symbol",
        "Gene_Type": "NCTC5055_Gene_Type",
        "Protein_Name": "NCTC5055_Protein_Name",
        "Protein_Accession": "NCTC5055_Protein_Accession",
        "Protein_Length": "NCTC5055_Protein_Length",
    }
)

merged = rbh.merge(
    pao1_keep,
    left_on="query_PAO1_to_NCTC5055",
    right_on="PAO1_Protein_Accession",
    how="left",
)

merged = merged.merge(
    nctc_keep,
    left_on="subject_PAO1_to_NCTC5055",
    right_on="NCTC5055_Protein_Accession",
    how="left",
)

merged.to_csv("PAO1_NCTC5055_ANNOTATED_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
merged.to_sql(
    "PAO1_NCTC5055_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("Saved CSV: PAO1_NCTC5055_ANNOTATED_LINKER.csv")
print("Imported into SQLite: PAO1_NCTC5055_ANNOTATED_LINKER")
print("Rows:", len(merged))
print("PAO1 annotation matches:", merged["PAO1_Protein_Name"].notna().sum())
print("NCTC5055 annotation matches:", merged["NCTC5055_Protein_Name"].notna().sum())