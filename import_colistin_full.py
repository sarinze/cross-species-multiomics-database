import pandas as pd
import sqlite3

file = "PA14_Colistin_Full_RNASeq_Data_Info.csv"

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
    "Colistin_RNAseq_Fold_Change",
    "Colistin_RNAseq_P_Value",
    "Colistin_RNAseq_BaseMean_Colistin",
    "Colistin_RNAseq_BaseMean_Control",
    "Colistin_Riboseq_Fold_Change",
    "Colistin_Riboseq_P_Value",
    "Colistin_Riboseq_BaseMean_Colistin",
    "Colistin_Riboseq_BaseMean_Control"
]

# remove empty rows
df = df.dropna(subset=["PA14_Original_Locus_Tag"])

print(df.head())
print("\nRows:", len(df))
print("\nColumns:")
print(df.columns.tolist())

conn = sqlite3.connect("pseudomonas_new.db")

df.to_sql(
    "PA14_COLISTIN_FULL_RNASEQ_INFO",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nImported clean full colistin RNA-seq/Ribo-seq table as: PA14_COLISTIN_FULL_RNASEQ_INFO")