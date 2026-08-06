import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

rbh = pd.read_sql_query("""
SELECT *
FROM KPPR1_MGH78578_RBH_LINKER
""", conn)

kppr1 = pd.read_sql_query("""
SELECT
    Locus_Tag AS KPPR1_Locus_Tag,
    Protein_Accession AS KPPR1_Protein_Accession,
    Protein_Name AS KPPR1_Protein_Name,
    Symbol AS KPPR1_Symbol
FROM KPPR1_Genome_Annotation
""", conn)

kppr1 = kppr1.sort_values("KPPR1_Locus_Tag")

kppr1 = kppr1.drop_duplicates(
    subset=["KPPR1_Protein_Accession"],
    keep="first"
)

mgh = pd.read_sql_query("""
SELECT
    Locus_Tag AS MGH78578_Locus_Tag,
    Protein_Accession AS MGH78578_Protein_Accession,
    Protein_Name AS MGH78578_Protein_Name,
    Symbol AS MGH78578_Symbol
FROM Kleb_MGH_78578_Genome_Annotation
WHERE Protein_Accession IS NOT NULL
""", conn)

mgh = mgh.sort_values("MGH78578_Locus_Tag")

mgh = mgh.drop_duplicates(
    subset=["MGH78578_Protein_Accession"],
    keep="first"
)

annotated = rbh.merge(
    kppr1,
    left_on="query_KPPR1_to_MGH78578",
    right_on="KPPR1_Protein_Accession",
    how="left"
)

annotated = annotated.merge(
    mgh,
    left_on="subject_KPPR1_to_MGH78578",
    right_on="MGH78578_Protein_Accession",
    how="left"
)

annotated.to_csv("KPPR1_MGH78578_ANNOTATED_LINKER.csv", index=False)

annotated.to_sql(
    "KPPR1_MGH78578_ANNOTATED_LINKER",
    conn,
    if_exists="replace",
    index=False
)

print("Annotated KPPR1-MGH78578 linker created")
print("Rows:", len(annotated))
print("Columns:", annotated.columns.tolist())

print("\nMissing KPPR1 locus tags:")
print(annotated["KPPR1_Locus_Tag"].isna().sum())

print("\nMissing MGH78578 locus tags:")
print(annotated["MGH78578_Locus_Tag"].isna().sum())

print(annotated.head())

conn.close()