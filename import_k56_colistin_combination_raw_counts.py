import pandas as pd
import sqlite3

# Load CSV file
df = pd.read_csv(
    "KPHS_K56_Colistin_Combination_Raw_Count_Table.csv"
)

# Show original structure
print("Original Columns:")
print(df.columns.tolist())
print("Column Count:", len(df.columns))
print("Rows:", len(df))

# Rename columns
df.columns = [
    "Gene_ID",
    "K56_Untreated",
    "K56_Colistin",
    "K56_Chemical3",
    "K56_CO_CHE3_Combination"
]

# Preview
print("\nFirst 5 Rows:")
print(df.head())

# Connect to SQLite database
conn = sqlite3.connect("pseudomonas_new.db")

# Import table
df.to_sql(
    "KPHS_K56_COLISTIN_COMBINATION_RAW_COUNTS",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nSUCCESS!")
print("Table Imported: KPHS_K56_COLISTIN_COMBINATION_RAW_COUNTS")
print("Rows:", len(df))
print("Columns:", len(df.columns))