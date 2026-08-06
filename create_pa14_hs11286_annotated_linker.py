import pandas as pd
import sqlite3

rbh = pd.read_csv("PA14_HS11286_RBH_Orthologs.csv")
pa14 = pd.read_csv("PA14_Genome_Annotation.csv")
hs = pd.read_csv("Kleb_HS11286_Genome_Annotation.csv")

# Keep useful PA14 columns
pa14_keep = pa14[
    [
        "Original_Locus_Tag",
        "New_Locus_Tag",
        "Gene_Name",
        "Symbol",
        "Gene_Type",
        "Protein_Name",
        "Protein_Accession",
        "Protein_Length",
    ]
].copy()

pa14_keep = pa14_keep.rename(
    columns={
        "Original_Locus_Tag": "PA14_Original_Locus_Tag",
        "New_Locus_Tag": "PA14_New_Locus_Tag",
        "Gene_Name": "PA14_Gene_Name",
        "Symbol": "PA14_Symbol",
        "Gene_Type": "PA14_Gene_Type",
        "Protein_Name": "PA14_Protein_Name",
        "Protein_Accession": "PA14_Protein_Accession",
        "Protein_Length": "PA14_Protein_Length",
    }
)

# Keep useful HS11286 columns
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

# Merge PA14 annotation
merged = rbh.merge(
    pa14_keep,
    left_on="query_PA14_to_HS11286",
    right_on="PA14_Protein_Accession",
    how="left",
)

# Merge HS11286 annotation
merged = merged.merge(
    hs_keep,
    left_on="subject_PA14_to_HS11286",
    right_on="HS11286_Protein_Accession",
    how="left",
)

merged.to_csv("PA14_HS11286_ANNOTATED_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
merged.to_sql(
    "PA14_HS11286_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("Saved CSV: PA14_HS11286_ANNOTATED_LINKER.csv")
print("Imported into SQLite: PA14_HS11286_ANNOTATED_LINKER")
print("Rows:", len(merged))
print("PA14 annotation matches:", merged["PA14_Protein_Name"].notna().sum())
print("HS11286 annotation matches:", merged["HS11286_Protein_Name"].notna().sum())