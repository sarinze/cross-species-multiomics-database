import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

# =========================
# 1. PA14 RNA-seq + Ribo-seq
# =========================

colistin = pd.read_sql_query("""
SELECT
    PA14_Original_Locus_Tag,
    Colistin_RNAseq_Fold_Change,
    Colistin_RNAseq_P_Value,
    Colistin_Riboseq_Fold_Change,
    Colistin_Riboseq_P_Value
FROM PA14_COLISTIN_FULL_RNASEQ_INFO
""", conn)

tobramycin = pd.read_sql_query("""
SELECT
    PA14_Original_Locus_Tag,
    Tobramycin_RNAseq_Fold_Change,
    Tobramycin_RNAseq_P_Value,
    Tobramycin_Riboseq_Fold_Change,
    Tobramycin_Riboseq_P_Value
FROM PA14_TOBRAMYCIN_FULL_RNASEQ_INFO
""", conn)

matrix = colistin.merge(
    tobramycin,
    on="PA14_Original_Locus_Tag",
    how="outer"
)

print("After RNA/Ribo merge:", matrix.shape)

# =========================
# 2. Deduplicate PA14-PAO1 linker BEFORE merging
# =========================

linker = pd.read_sql_query("""
SELECT
    PA14_Original_Locus_Tag,
    PAO1_Locus_Tag,
    PAO1_Gene_Name,
    PAO1_Symbol,
    BITSCORE,
    PERCENT_IDENTITY
FROM PA14_PAO1_CLEAN_LINKER
""", conn)

print("Linker before deduplication:", linker.shape)

linker = linker.sort_values(
    by=["PA14_Original_Locus_Tag", "BITSCORE"],
    ascending=[True, False]
)

linker = linker.drop_duplicates(
    subset=["PA14_Original_Locus_Tag"],
    keep="first"
)

print("Linker after deduplication:", linker.shape)

matrix = matrix.merge(
    linker,
    on="PA14_Original_Locus_Tag",
    how="left"
)

print("After adding PAO1 linker:", matrix.shape)

# =========================
# 3. PAO1 Tobramycin Tn-seq
# =========================

pao1_tobra_tnseq = pd.read_sql_query("""
SELECT
    PAO1_Locus_Tag,
    Selection_Ratio AS PAO1_Tobramycin_TnSeq_Selection_Ratio,
    Pre_Growth_Hits AS PAO1_Tobramycin_Pre_Growth_Hits,
    Growth_With_Tobramycin_Hits AS PAO1_Tobramycin_Growth_With_Tobramycin_Hits
FROM PAO1_Tobramycin_TnSeq
""", conn)

pao1_tobra_tnseq = pao1_tobra_tnseq.drop_duplicates(
    subset=["PAO1_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pao1_tobra_tnseq,
    on="PAO1_Locus_Tag",
    how="left"
)

print("After adding Tobramycin Tn-seq:", matrix.shape)

# =========================
# 4. PAO1 Persister Tn-seq
# =========================

pao1_persister = pd.read_sql_query("""
SELECT
    PAO1_Locus_Tag,
    Mean_Survival_Index AS PAO1_Persister_TnSeq_Mean_Survival_Index
FROM PAO1_Persister_TnSeq_All_SI
""", conn)

pao1_persister = pao1_persister.drop_duplicates(
    subset=["PAO1_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pao1_persister,
    on="PAO1_Locus_Tag",
    how="left"
)

print("After adding Persister Tn-seq:", matrix.shape)

# =========================
# 5. Save final matrix
# =========================

matrix.to_csv("PA14_RNA_RIBO_PAO1_TNSEQ_MATRIX.csv", index=False)

print("\nMatrix created successfully")
print("Final rows:", len(matrix))
print("Final columns:")
print(matrix.columns.tolist())
print(matrix.head())

print("\nMissing PAO1 Tn-seq values:")
print(matrix[
    [
        "PAO1_Tobramycin_TnSeq_Selection_Ratio",
        "PAO1_Persister_TnSeq_Mean_Survival_Index"
    ]
].isna().sum())

conn.close()