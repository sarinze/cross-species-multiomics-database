import pandas as pd
import sqlite3

rbh = pd.read_csv("PAO1_HS11286_RBH_Orthologs.csv")
pao1 = pd.read_csv("PAO1_Genome_Annotation.csv")
hs = pd.read_csv("Kleb_HS11286_Genome_Annotation.csv")

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

hs_keep = hs[
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

hs_keep = hs_keep.rename(
    columns={
        "Locus_Tag": "HS11286_Locus_Tag",
        "Gene_Name": "HS11286_Gene_Name",
        "Symbol": "HS11286_Symbol",
        "Gene_Type": "HS11286_Gene_Type",
        "Protein_Name": "HS11286_Protein_Name",
        "Protein_Accession": "HS11286_Protein_Accession",
        "Protein_Length": "HS11286_Protein_Length",
    }
)

merged = rbh.merge(
    pao1_keep,
    left_on="query_PAO1_to_HS11286",
    right_on="PAO1_Protein_Accession",
    how="left",
)

merged = merged.merge(
    hs_keep,
    left_on="subject_PAO1_to_HS11286",
    right_on="HS11286_Protein_Accession",
    how="left",
)

merged.to_csv("PAO1_HS11286_ANNOTATED_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
merged.to_sql(
    "PAO1_HS11286_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("Saved CSV: PAO1_HS11286_ANNOTATED_LINKER.csv")
print("Imported into SQLite: PAO1_HS11286_ANNOTATED_LINKER")
print("Rows:", len(merged))
print("PAO1 annotation matches:", merged["PAO1_Protein_Name"].notna().sum())
print("HS11286 annotation matches:", merged["HS11286_Protein_Name"].notna().sum())