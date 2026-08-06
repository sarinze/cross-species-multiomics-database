import sqlite3
import pandas as pd

conn = sqlite3.connect("pseudomonas_new.db")

tables_to_check = [
    "PAO1_Tobramycin_TnSeq",
    "PAO1_Persister_TnSeq_All_SI",
    "PAO1_Persister_TnSeq_10fold_Decrease_SI",
    "PAO1_Persister_TnSeq_10fold_Increase_SI",
    "PA14_PAO1_CLEAN_LINKER"
]

for table in tables_to_check:
    print("\n==============================")
    print(table)
    print("==============================")

    df = pd.read_sql_query(f"SELECT * FROM {table} LIMIT 5", conn)

    print("Columns:")
    print(df.columns.tolist())

    print("\nFirst 5 rows:")
    print(df)

conn.close()