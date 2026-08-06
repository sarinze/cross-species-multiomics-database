import pandas as pd
import sqlite3

rbh = pd.read_csv("PA14_MGH78578_RBH_Orthologs.csv")
pa14 = pd.read_csv("PA14_Genome_Annotation.csv")
mgh = pd.read_csv("Kleb_MGH_78578_Genome_Annotation.csv")

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

# Keep useful MGH78578 columns
mgh_keep = mgh[
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

mgh_keep = mgh_keep.rename(
    columns={
        "Locus_Tag": "MGH78578_Locus_Tag",
        "Gene_Name": "MGH78578_Gene_Name",
        "Symbol": "MGH78578_Symbol",
        "Gene_Type": "MGH78578_Gene_Type",
        "Protein_Name": "MGH78578_Protein_Name",
        "Protein_Accession": "MGH78578_Protein_Accession",
        "Protein_Length": "MGH78578_Protein_Length",
    }
)

# Merge PA14 annotation
merged = rbh.merge(
    pa14_keep,
    left_on="queryPA14_to_MGH78578",
    right_on="PA14_Protein_Accession",
    how="left",
)

# Merge MGH78578 annotation
merged = merged.merge(
    mgh_keep,
    left_on="subjectPA14_to_MGH78578",
    right_on="MGH78578_Protein_Accession",
    how="left",
)

merged.to_csv("PA14_MGH78578_ANNOTATED_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
merged.to_sql(
    "PA14_MGH78578_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("Saved CSV: PA14_MGH78578_ANNOTATED_LINKER.csv")
print("Imported into SQLite: PA14_MGH78578_ANNOTATED_LINKER")
print("Rows:", len(merged))
print("PA14 annotation matches:", merged["PA14_Protein_Name"].notna().sum())
print("MGH78578 annotation matches:", merged["MGH78578_Protein_Name"].notna().sum())