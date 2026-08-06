import pandas as pd
import sqlite3

# ------------------------------------------------------------
# PA14 vs NCTC5055 reciprocal best hit linker script
# Purpose:
# This script identifies reciprocal best-hit protein matches
# between Pseudomonas aeruginosa PA14 and
# Klebsiella pneumoniae NCTC5055.
# ------------------------------------------------------------

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

# Load BLASTP output files
pa14_to_nctc = pd.read_csv(
    "PA14_vs_NCTC5055_REFSEQ_CLEAN.tsv",
    sep="\t",
    names=cols,
)

nctc_to_pa14 = pd.read_csv(
    "NCTC5055_vs_PA14_REFSEQ_CLEAN.tsv",
    sep="\t",
    names=cols,
)

# Select the best hit for each PA14 protein
pa14_best = (
    pa14_to_nctc.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

# Select the best hit for each NCTC5055 protein
nctc_best = (
    nctc_to_pa14.sort_values("bitscore", ascending=False)
    .drop_duplicates("query")
)

# Find reciprocal best hits
rbh = pa14_best.merge(
    nctc_best,
    left_on=["query", "subject"],
    right_on=["subject", "query"],
    suffixes=("_PA14_to_NCTC5055", "_NCTC5055_to_PA14"),
)

# Save output as CSV
rbh.to_csv("PA14_NCTC5055_RBH_Orthologs.csv", index=False)

# Save output into SQLite database
conn = sqlite3.connect("pseudomonas_new.db")

rbh.to_sql(
    "PA14_NCTC5055_RBH_LINKER",
    conn,
    if_exists="replace",
    index=False,
)

conn.close()

# Print summary
print("Reciprocal Best Hits:", len(rbh))
print("Saved CSV: PA14_NCTC5055_RBH_Orthologs.csv")
print("Imported into SQLite: PA14_NCTC5055_RBH_LINKER")