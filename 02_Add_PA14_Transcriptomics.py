import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"
INPUT_FILE = "Master_Ortholog_Backbone.csv"
OUTPUT_FILE = "Master_Transcriptomics_Matrix.csv"

conn = sqlite3.connect(DB_NAME)

matrix = pd.read_csv(INPUT_FILE)

# =========================
# Add PA14 Colistin RNA-seq/Ribo-seq
# =========================

pa14_colistin = pd.read_sql_query("""
SELECT
    PA14_Original_Locus_Tag,
    Colistin_RNAseq_Fold_Change,
    Colistin_RNAseq_P_Value,
    Colistin_RNAseq_BaseMean_Colistin,
    Colistin_RNAseq_BaseMean_Control,
    Colistin_Riboseq_Fold_Change,
    Colistin_Riboseq_P_Value,
    Colistin_Riboseq_BaseMean_Colistin,
    Colistin_Riboseq_BaseMean_Control
FROM PA14_COLISTIN_FULL_RNASEQ_INFO
""", conn)

pa14_colistin = pa14_colistin.drop_duplicates(
    subset=["PA14_Original_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pa14_colistin,
    left_on="Pseudomonas_Original_Locus_Tag",
    right_on="PA14_Original_Locus_Tag",
    how="left"
)

matrix = matrix.drop(columns=["PA14_Original_Locus_Tag"])

print("After adding PA14 colistin:", matrix.shape)

# =========================
# Add PA14 Tobramycin RNA-seq/Ribo-seq
# =========================

pa14_tobramycin = pd.read_sql_query("""
SELECT
    PA14_Original_Locus_Tag,
    Tobramycin_RNAseq_Fold_Change,
    Tobramycin_RNAseq_P_Value,
    Tobramycin_RNAseq_BaseMean_Tobramycin,
    Tobramycin_RNAseq_BaseMean_Control,
    Tobramycin_Riboseq_Fold_Change,
    Tobramycin_Riboseq_P_Value,
    Tobramycin_Riboseq_BaseMean_Tobramycin,
    Tobramycin_Riboseq_BaseMean_Control
FROM PA14_TOBRAMYCIN_FULL_RNASEQ_INFO
""", conn)

pa14_tobramycin = pa14_tobramycin.drop_duplicates(
    subset=["PA14_Original_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pa14_tobramycin,
    left_on="Pseudomonas_Original_Locus_Tag",
    right_on="PA14_Original_Locus_Tag",
    how="left"
)

matrix = matrix.drop(columns=["PA14_Original_Locus_Tag"])

print("After adding PA14 tobramycin:", matrix.shape)

matrix.to_csv(OUTPUT_FILE, index=False)

print("\nPA14 transcriptomics added successfully")
print("Output file:", OUTPUT_FILE)
print("Rows:", len(matrix))
print("Columns:", len(matrix.columns))

conn.close()