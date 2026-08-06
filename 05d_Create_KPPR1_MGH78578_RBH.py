import pandas as pd
import sqlite3

cols = [
    "query", "subject", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

forward = pd.read_csv("KPPR1_vs_MGH78578.tsv", sep="\t", names=cols)
reverse = pd.read_csv("MGH78578_vs_KPPR1.tsv", sep="\t", names=cols)

kppr1_best = (
    forward.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

mgh_best = (
    reverse.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

rbh = kppr1_best.merge(
    mgh_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("_KPPR1_to_MGH78578", "_MGH78578_to_KPPR1")
)

rbh.to_csv("KPPR1_MGH78578_RBH_Orthologs.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
rbh.to_sql(
    "KPPR1_MGH78578_RBH_LINKER",
    conn,
    if_exists="replace",
    index=False
)
conn.close()

print("KPPR1-MGH78578 RBH linker created")
print("Rows:", len(rbh))
print("Columns:", rbh.columns.tolist())
print(rbh.head())