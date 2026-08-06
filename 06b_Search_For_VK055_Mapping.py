import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"

conn = sqlite3.connect(DB_NAME)

tables = pd.read_sql_query(
    "SELECT name FROM sqlite_master WHERE type='table';",
    conn
)["name"].tolist()

for table in tables:
    try:
        df = pd.read_sql_query(f"SELECT * FROM {table} LIMIT 5", conn)
        cols = df.columns.tolist()

        # Search column names
        matching_cols = [c for c in cols if "VK055" in c.upper() or "GENE" in c.upper() or "LOCUS" in c.upper()]

        if matching_cols:
            print("\n==============================")
            print(table)
            print("==============================")
            print("Relevant columns:", matching_cols)
            print(df[matching_cols].head())

    except Exception as e:
        pass

conn.close()