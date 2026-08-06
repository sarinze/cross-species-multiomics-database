import pandas as pd
import sqlite3

cols = [
    "query", "subject", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

pao1_to_nctc = pd.read_csv(
    "PAO1_vs_NCTC5055_CLEAN.tsv",
    sep="\t",
    names=cols,
)

nctc_to_pao1 = pd.read_csv(
    "NCTC5055_vs_PAO1_CLEAN.tsv",
    sep="\t",
    names=cols,
)

pao1_best = (
    pao1_to_nctc.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

nctc_best = (
    nctc_to_pao1.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

rbh = pao1_best.merge(
    nctc_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("_PAO1_to_NCTC5055", "_NCTC5055_to_PAO1"),
)

rbh.to_csv("PAO1_NCTC5055_RBH_Orthologs.csv", index=False)

conn = sqlite3.connect("pseudomonas_new.db")
rbh.to_sql("PAO1_NCTC5055_RBH_LINKER", conn, if_exists="replace", index=False)
conn.close()

print("Reciprocal Best Hits:", len(rbh))
print("Saved CSV: PAO1_NCTC5055_RBH_Orthologs.csv")
print("Imported into SQLite: PAO1_NCTC5055_RBH_LINKER")