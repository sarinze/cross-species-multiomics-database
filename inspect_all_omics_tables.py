import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

tables_to_check = [
    "PA14_COLISTIN_RNASEQ",
    "PA14_TOBRAMYCIN_RNASEQ",
    "PA14_COLISTIN_FULL_RNASEQ_INFO",
    "PA14_TOBRAMYCIN_FULL_RNASEQ_INFO",

    "PAO1_Tobramycin_TnSeq",
    "PAO1_Persister_TnSeq_All_SI",
    "PAO1_Persister_TnSeq_10fold_Decrease_SI",
    "PAO1_Persister_TnSeq_10fold_Increase_SI",

    "PAO1_PolymycinB_Rifampicin_Annotation",
    "PAO1_PolymycinB_Rifampicin_Content",

    "Klebsiella_TnSeq",
    "KPHS_K56_COLISTIN_COMBINATION_DEG",
    "KPHS_K56_COLISTIN_COMBINATION_RAW_COUNTS"
]

for table in tables_to_check:
    print("\n" + "=" * 80)
    print(table)
    print("=" * 80)

    try:
        df = pd.read_sql_query(f"SELECT * FROM {table} LIMIT 5", conn)

        print("Columns:")
        print(df.columns.tolist())

        print("\nShape:")
        count = pd.read_sql_query(f"SELECT COUNT(*) AS rows FROM {table}", conn)
        print(count)

        print("\nFirst 5 rows:")
        print(df)

    except Exception as e:
        print("Could not inspect table:")
        print(e)

conn.close()