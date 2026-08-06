import pandas as pd
import sqlite3
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
    "bitscore"
]
forward = pd.read_csv("PA14_vs_MGH78578_REFSEQ_CLEAN.tsv", sep="\t", names=cols)
reverse = pd.read_csv("MGH78578_vs_PA14_REFSEQ_CLEAN.tsv", sep="\t", names=cols)

forward_best = (
forward.sort_values("bitscore", ascending=False)
.drop_duplicates("query")
)
reverse_best = (
reverse.sort_values("bitscore", ascending=False)
.drop_duplicates("query")
)
rbh = forward_best.merge(
    reverse_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("PA14_to_MGH78578", "MGH78578_to_PA14")
)
print("Reciprocal Best Hits:", len(rbh))

rbh.to_csv("PA14_MGH78578_RBH_Orthologs.csv", index=False)
conn = sqlite3.connect("pseudomonas_new.db")

rbh.to_sql("PA14_MGH78578_RBH_LINKER", conn, if_exists="replace", index=False)

conn.close()

print("Saved CSV: PA14_MGH78578_RBH_Orthologs.csv")
print("Imported into SQLite: PA14_MGH78578_RBH_LINKER")