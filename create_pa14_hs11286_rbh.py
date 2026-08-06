import pandas as pd
import sqlite3

cols = [
    "query", "subject", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

pa14_to_hs11286 = pd.read_csv(
    "PA14_vs_HS11286_REFSEQ_CLEAN.tsv",
    sep="\t",
    names=cols,
)

hs11286_to_pa14 = pd.read_csv(
    "HS11286_vs_PA14_REFSEQ_CLEAN.tsv",
    sep="\t",
    names=cols,
)

pa14_best = (
    pa14_to_hs11286.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

hs11286_best = (
    hs11286_to_pa14.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

rbh = pa14_best.merge(
    hs11286_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("_PA14_to_HS11286", "_HS11286_to_PA14"),
)

rbh.to_csv("PA14_HS11286_RBH_Orthologs.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
rbh.to_sql("PA14_HS11286_RBH_LINKER", conn, if_exists="replace", index=False)
conn.close()

print("Reciprocal Best Hits:", len(rbh))
print("Saved CSV: PA14_HS11286_RBH_Orthologs.csv")
print("Imported into SQLite: PA14_HS11286_RBH_LINKER")