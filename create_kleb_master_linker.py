import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

hs_nctc = pd.read_sql("SELECT * FROM Kleb_HS11286_NCTC5055_PROT_RBH_LINKER", conn)
hs_mgh = pd.read_sql("SELECT * FROM Kleb_HS11286_MGH78578_PROT_RBH_LINKER", conn)

hs_annot = pd.read_sql("SELECT * FROM Kleb_HS11286_Genome_Annotation", conn)
nctc_annot = pd.read_sql("SELECT * FROM NCTC5055_Genome_Annotation", conn)
mgh_annot = pd.read_sql("SELECT * FROM Kleb_MGH_78578_Genome_Annotation", conn)

# Remove blank protein accessions and duplicate protein accessions before merging
hs_annot = hs_annot.dropna(subset=["Protein_Accession"]).drop_duplicates("Protein_Accession")
nctc_annot = nctc_annot.dropna(subset=["Protein_Accession"]).drop_duplicates("Protein_Accession")
mgh_annot = mgh_annot.dropna(subset=["Protein_Accession"]).drop_duplicates("Protein_Accession")

# Rename annotation columns so each strain is clearly labelled
hs_annot = hs_annot.add_prefix("HS11286_")
nctc_annot = nctc_annot.add_prefix("NCTC5055_")
mgh_annot = mgh_annot.add_prefix("MGH78578_")

master = hs_nctc.merge(
    hs_mgh,
    on="Kleb_HS11286_SEQUENCE",
    how="outer",
    suffixes=("_HS11286_NCTC5055", "_HS11286_MGH78578")
)

master = master.merge(
    hs_annot,
    left_on="Kleb_HS11286_SEQUENCE",
    right_on="HS11286_Protein_Accession",
    how="left"
)

master = master.merge(
    nctc_annot,
    left_on="NCTC5055_SEQUENCE",
    right_on="NCTC5055_Protein_Accession",
    how="left"
)

master = master.merge(
    mgh_annot,
    left_on="MGH78578_SEQUENCE",
    right_on="MGH78578_Protein_Accession",
    how="left"
)

master.to_csv("KLEBSIELLA_MASTER_LINKER.csv", index=False)

master.to_sql(
    "KLEBSIELLA_MASTER_LINKER",
    conn,
    if_exists="replace",
    index=False
)

print("Saved CSV: KLEBSIELLA_MASTER_LINKER.csv")
print("Imported into SQLite as: KLEBSIELLA_MASTER_LINKER")
print("Rows:", len(master))
print("Columns:")
print(master.columns.tolist())

conn.close()