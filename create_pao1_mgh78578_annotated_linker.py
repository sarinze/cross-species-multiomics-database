import pandas as pd
import sqlite3

rbh = pd.read_csv("PAO1_MGH78578_RBH_Orthologs.csv")
pao1 = pd.read_csv("PAO1_Genome_Annotation.csv")
mgh = pd.read_csv("Kleb_MGH_78578_Genome_Annotation.csv")

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

merged = rbh.merge(
    pao1_keep,
    left_on="query_PAO1_to_MGH78578",
    right_on="PAO1_Protein_Accession",
    how="left",
)

merged = merged.merge(
    mgh_keep,
    left_on="subject_PAO1_to_MGH78578",
    right_on="MGH78578_Protein_Accession",
    how="left",
)

merged.to_csv("PAO1_MGH78578_ANNOTATED_LINKER.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
merged.to_sql(
    "PAO1_MGH78578_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False,
)
conn.close()

print("Saved CSV: PAO1_MGH78578_ANNOTATED_LINKER.csv")
print("Imported into SQLite: PAO1_MGH78578_ANNOTATED_LINKER")
print("Rows:", len(merged))
print("PAO1 annotation matches:", merged["PAO1_Protein_Name"].notna().sum())
print("MGH78578 annotation matches:", merged["MGH78578_Protein_Name"].notna().sum())