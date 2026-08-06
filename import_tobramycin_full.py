import pandas as pd
import sqlite3

file = "PA14_Tobramycin_Full_RNASeq_Data_Info.csv"

df = pd.read_csv(file, header=3)

# keep only the real gene-level data columns
df = df.iloc[:, 0:15]

df.columns = [
    "Gene",
    "PA14_Original_Locus_Tag",
    "Start",
    "End",
    "Strand",
    "PseudoCAP_Category",
    "Product_Name",
    "Tobramycin_RNAseq_Fold_Change",
    "Tobramycin_RNAseq_P_Value",
    "Tobramycin_RNAseq_BaseMean_Tobramycin",
    "Tobramycin_RNAseq_BaseMean_Control",
    "Tobramycin_Riboseq_Fold_Change",
    "Tobramycin_Riboseq_P_Value",
    "Tobramycin_Riboseq_BaseMean_Tobramycin",
    "Tobramycin_Riboseq_BaseMean_Control"
]

# remove empty rows
df = df.dropna(subset=["PA14_Original_Locus_Tag"])

print(df.head())
print("\nRows:", len(df))
print("\nColumns:")
print(df.columns.tolist())

conn = sqlite3.connect("pseudomonas_new.db")

df.to_sql(
    "PA14_TOBRAMYCIN_FULL_RNASEQ_INFO",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nImported clean full tobramycin RNA-seq/Ribo-seq table as: PA14_TOBRAMYCIN_FULL_RNASEQ_INFO")