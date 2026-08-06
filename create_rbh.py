import pandas as pd
cols = [
    "query",
    "subject",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]
PAO1_PA14 = pd.read_csv(
    "PAO1_vs_PA14.tsv",
    sep="\t",
    names=cols
)
PA14_PAO1 = pd.read_csv(
    "PA14_vs_PAO1.tsv",
    sep="\t",
    names=cols
)
print("PAO1 vs PA14:")
print(PAO1_PA14.head())
print("\nPA14 vs PAO1:")
print(PA14_PAO1.head())
# Select best hits for each query based on highest bitscore
PAO1_best = (
    PAO1_PA14.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)
PA14_best = (
    PA14_PAO1.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)
# Find reciprocal best hits
rbh = PAO1_best.merge(
    PA14_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("_PAO1_to_PA14", "_PA14_to_PAO1")
)
print(rbh.columns.tolist())
orthologs = rbh[[
    "query_PAO1_to_PA14",
    "subject_PAO1_to_PA14",
    "pident_PAO1_to_PA14",
    "bitscore_PAO1_to_PA14"
]]
orthologs.columns = [
    "PAO1_PROTEIN",
    "PA14_PROTEIN",
    "PERCENT_IDENTITY",
    "BITSCORE"
]
print("Number of reciprocal best-hit orthologs:",len(orthologs))
print(orthologs.head())
orthologs.to_csv("PAO1_PA14_RBH_Orthologs.csv", index=False)
print("Saved file: PAO1_PA14_RBH_Orthologs.csv")


pao1_annot = pd.read_csv("PAO1_Genome_Annotation.csv")
pa14_annot = pd.read_csv("PA14_Genome_Annotation.csv")

final_linker = orthologs.merge(
    pao1_annot,
    left_on="PAO1_PROTEIN",
    right_on="Protein_Accession",
    how="left",
    suffixes=("_PAO1", "_PA14")
    )
final_linker = final_linker.merge(
        pa14_annot,
        left_on="PA14_PROTEIN",
        right_on="Protein_Accession",
        how="left",
        suffixes=("_PAO1", "_PA14")
    )
final_linker.to_csv("PAO1_PA14_Final_Locus_Linker.csv", index=False)
print("Saved final linker table: PAO1_PA14_Final_Locus_Linker.csv")
print(final_linker.head())


import sqlite3
conn = sqlite3.connect("Pseudomonas_new.db")
final_linker.to_sql("PAO1_PA14_RBH_LINKER", conn, if_exists="replace", index=False)
conn.close()
print("imported into SQLite as table: PAO1_PA14_RBH_LINKER")