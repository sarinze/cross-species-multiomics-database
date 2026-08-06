import pandas as pd
import sqlite3

file = "PA14_Tobramycin_RNASeq_Data.csv"

df = pd.read_csv(file, header=None, skiprows=4)

df = df.iloc[:, 0:3]

df.columns = [
    "PA14_Original_Locus_Tag",
    "Tobramycin_RNAseq_Fold_Change",
    "Tobramycin_Riboseq_Fold_Change"
]

df = df.dropna(subset=["PA14_Original_Locus_Tag"])

print(df.head())
print("\nRows:", len(df))
print("\nColumns:", df.columns.tolist())

conn = sqlite3.connect("pseudomonas_new.db")

df.to_sql(
    "PA14_TOBRAMYCIN_RNASEQ",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nImported clean RNA-seq table as: PA14_TOBRAMYCIN_RNASEQ")