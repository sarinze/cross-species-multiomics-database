import pandas as pd
import sqlite3

file = "PA14_Colistin_RNASeq_Data.csv"

df = pd.read_csv(file, header=None, skiprows=4)

# keep only the first 3 useful columns
df = df.iloc[:, 0:3]

# rename columns clearly
df.columns = [
    "PA14_Original_Locus_Tag",
    "Colistin_RNAseq_Fold_Change",
    "Colistin_Riboseq_Fold_Change"
]

# remove empty rows
df = df.dropna(subset=["PA14_Original_Locus_Tag"])

print(df.head())
print("\nRows:", len(df))
print("\nColumns:", df.columns.tolist())

conn = sqlite3.connect("pseudomonas_new.db")

df.to_sql(
    "PA14_COLISTIN_RNASEQ",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nImported clean RNA-seq table as: PA14_COLISTIN_RNASEQ")