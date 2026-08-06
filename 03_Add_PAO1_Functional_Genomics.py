import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"
INPUT_FILE = "Master_Transcriptomics_Matrix.csv"
OUTPUT_FILE = "Master_Functional_Genomics_Matrix.csv"

conn = sqlite3.connect(DB_NAME)

matrix = pd.read_csv(INPUT_FILE)

# =========================
# PAO1 Tobramycin Tn-seq
# =========================

pao1_tobra = pd.read_sql_query("""
SELECT
    PAO1_Locus_Tag,
    Pre_Growth_Hits AS PAO1_Tobramycin_Pre_Growth_Hits,
    Pre_Growth_Reads AS PAO1_Tobramycin_Pre_Growth_Reads,
    Growth_No_Tobramycin_Hits AS PAO1_Tobramycin_No_Tobra_Hits,
    Growth_No_Tobramycin_Reads AS PAO1_Tobramycin_No_Tobra_Reads,
    Growth_With_Tobramycin_Hits AS PAO1_Tobramycin_With_Tobra_Hits,
    Growth_With_Tobramycin_Reads AS PAO1_Tobramycin_With_Tobra_Reads,
    Selection_Ratio AS PAO1_Tobramycin_Selection_Ratio
FROM PAO1_Tobramycin_TnSeq
""", conn)

pao1_tobra = pao1_tobra.drop_duplicates(
    subset=["PAO1_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pao1_tobra,
    left_on="Pseudomonas_Locus_Tag",
    right_on="PAO1_Locus_Tag",
    how="left"
)

matrix = matrix.drop(columns=["PAO1_Locus_Tag"])

print("After adding PAO1 Tobramycin Tn-seq:", matrix.shape)

# =========================
# PAO1 Persister All SI
# =========================

pao1_persister_all = pd.read_sql_query("""
SELECT
    PAO1_Locus_Tag,
    Read_Count AS PAO1_Persister_All_Read_Count,
    Rep1_Survival_Index AS PAO1_Persister_All_Rep1_SI,
    Rep2_Survival_Index AS PAO1_Persister_All_Rep2_SI,
    Rep3_Survival_Index AS PAO1_Persister_All_Rep3_SI,
    Mean_Survival_Index AS PAO1_Persister_All_Mean_SI
FROM PAO1_Persister_TnSeq_All_SI
""", conn)

pao1_persister_all = pao1_persister_all.drop_duplicates(
    subset=["PAO1_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pao1_persister_all,
    left_on="Pseudomonas_Locus_Tag",
    right_on="PAO1_Locus_Tag",
    how="left"
)

matrix = matrix.drop(columns=["PAO1_Locus_Tag"])

print("After adding PAO1 Persister All SI:", matrix.shape)

# =========================
# PAO1 Persister 10fold Decrease SI
# =========================

pao1_persister_decrease = pd.read_sql_query("""
SELECT
    PAO1_Locus_Tag,
    Mean_Survival_Index AS PAO1_Persister_10fold_Decrease_Mean_SI
FROM PAO1_Persister_TnSeq_10fold_Decrease_SI
""", conn)

pao1_persister_decrease = pao1_persister_decrease.drop_duplicates(
    subset=["PAO1_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pao1_persister_decrease,
    left_on="Pseudomonas_Locus_Tag",
    right_on="PAO1_Locus_Tag",
    how="left"
)

matrix = matrix.drop(columns=["PAO1_Locus_Tag"])

print("After adding PAO1 Persister 10fold Decrease SI:", matrix.shape)

# =========================
# PAO1 Persister 10fold Increase SI
# =========================

pao1_persister_increase = pd.read_sql_query("""
SELECT
    PAO1_Locus_Tag,
    Mean_Survival_Index AS PAO1_Persister_10fold_Increase_Mean_SI
FROM PAO1_Persister_TnSeq_10fold_Increase_SI
""", conn)

pao1_persister_increase = pao1_persister_increase.drop_duplicates(
    subset=["PAO1_Locus_Tag"],
    keep="first"
)

matrix = matrix.merge(
    pao1_persister_increase,
    left_on="Pseudomonas_Locus_Tag",
    right_on="PAO1_Locus_Tag",
    how="left"
)

matrix = matrix.drop(columns=["PAO1_Locus_Tag"])

print("After adding PAO1 Persister 10fold Increase SI:", matrix.shape)

matrix.to_csv(OUTPUT_FILE, index=False)

print("\nPAO1 functional genomics added successfully")
print("Output file:", OUTPUT_FILE)
print("Rows:", len(matrix))
print("Columns:", len(matrix.columns))

conn.close()