import pandas as pd
import sqlite3

rbh = pd.read_csv("PA14_NCTC5055_RBH_Orthologs.csv")
pa14 = pd.read_csv("PA14_Genome_Annotation.csv")
nctc = pd.read_csv("NCTC5055_Genome_Annotation.csv")

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
    pa14_keep,
    left_on="query_PA14_to_NCTC5055",
    right_on="PA14_Protein_Accession",
    how="left",
)

merged = merged.merge(
    nctc_keep,
    left_on="subject_PA14_to_NCTC5055",
    right_on="NCTC5055_Protein_Accession",
    how="left",
)

merged.to_csv("PA14_NCTC5055_ANNOTATED_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
merged.to_sql(
    "PA14_NCTC5055_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("Saved CSV: PA14_NCTC5055_ANNOTATED_LINKER.csv")
print("Imported into SQLite: PA14_NCTC5055_ANNOTATED_LINKER")
print("Rows:", len(merged))
print("PA14 annotation matches:", merged["PA14_Protein_Name"].notna().sum())
print("NCTC5055 annotation matches:", merged["NCTC5055_Protein_Name"].notna().sum())