import pandas as pd
import sqlite3

# Load Excel file
df = pd.read_excel(
    "KPHS_K56_Colistin_Combination_DEG_Table.xlsx",
    header=4
)

# Check original structure
print("Original Columns:", len(df.columns))
print("Original Rows:", len(df))

# Rename columns
df.columns = [

    "Gene",

    # K56 vs Colistin
    "K56_vs_Colistin_Mean_Untreated",
    "K56_vs_Colistin_Mean_Colistin",
    "K56_vs_Colistin_Log2FC",
    "K56_vs_Colistin_Pooled_SD",
    "K56_vs_Colistin_Z_Statistic",
    "K56_vs_Colistin_P_Value",
    "K56_vs_Colistin_Q_Value",

    # K56 vs Chemical3
    "K56_vs_Chemical3_Mean_Untreated",
    "K56_vs_Chemical3_Mean_Chemical3",
    "K56_vs_Chemical3_Log2FC",
    "K56_vs_Chemical3_Pooled_SD",
    "K56_vs_Chemical3_Z_Statistic",
    "K56_vs_Chemical3_P_Value",
    "K56_vs_Chemical3_Q_Value",

    # K56 vs Combination
    "K56_vs_Combination_Mean_Untreated",
    "K56_vs_Combination_Mean_Combination",
    "K56_vs_Combination_Log2FC",
    "K56_vs_Combination_Pooled_SD",
    "K56_vs_Combination_Z_Statistic",
    "K56_vs_Combination_P_Value",
    "K56_vs_Combination_Q_Value",

    # Colistin vs Chemical3
    "Colistin_vs_Chemical3_Mean_Colistin",
    "Colistin_vs_Chemical3_Mean_Chemical3",
    "Colistin_vs_Chemical3_Log2FC",
    "Colistin_vs_Chemical3_Pooled_SD",
    "Colistin_vs_Chemical3_Z_Statistic",
    "Colistin_vs_Chemical3_P_Value",
    "Colistin_vs_Chemical3_Q_Value",

    # Colistin vs Combination
    "Colistin_vs_Combination_Mean_Colistin",
    "Colistin_vs_Combination_Mean_Combination",
    "Colistin_vs_Combination_Log2FC",
    "Colistin_vs_Combination_Pooled_SD",
    "Colistin_vs_Combination_Z_Statistic",
    "Colistin_vs_Combination_P_Value",
    "Colistin_vs_Combination_Q_Value",

    # Chemical3 vs Combination
    "Chemical3_vs_Combination_Mean_Chemical3",
    "Chemical3_vs_Combination_Mean_Combination",
    "Chemical3_vs_Combination_Log2FC",
    "Chemical3_vs_Combination_Pooled_SD",
    "Chemical3_vs_Combination_Z_Statistic",
    "Chemical3_vs_Combination_P_Value",
    "Chemical3_vs_Combination_Q_Value"
]

# Display validation information
print("\nRenamed Columns:")
for col in df.columns:
    print(col)

print("\nFinal Column Count:", len(df.columns))
print("Final Row Count:", len(df))

print("\nFirst 5 Rows:")
print(df.head())

# Connect to SQLite database
conn = sqlite3.connect("pseudomonas_new.db")

# Import table
df.to_sql(
    "KPHS_K56_COLISTIN_COMBINATION_DEG",
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\nSUCCESS!")
print("Table Imported: KPHS_K56_COLISTIN_COMBINATION_DEG")
print("Rows:", len(df))
print("Columns:", len(df.columns))