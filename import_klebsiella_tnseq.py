import pandas as pd
import sqlite3

# ============================
# Files
# ============================

csv_file = "Klebsiella_TnSeq_Table.csv"
db_file = "pseudomonas_new.db"
table_name = "Klebsiella_TnSeq"

# ============================
# Read CSV
# ============================

df = pd.read_csv(csv_file, low_memory=False)

# ============================
# Standardize column names
# ============================

df.columns = [
    "Gene_ID",
    "Gene_Name",
    "Gene_Length",
    "Total_Insertions",
    "Unique_Insertions",
    "Mean_Input",
    "Mean_Output",
    "Fold_Change",
    "Log2_Fold_Change_Output_Input",
    "P_Value",
    "Adjusted_P_Value",
    "Primary_KEGG_Annotation",
    "Secondary_KEGG_Annotation"
]

# ============================
# Quality Control
# ============================

print("\n========== QUALITY CONTROL ==========")

print("Rows:", len(df))
print("Columns:", len(df.columns))

print("\nMissing Gene_ID values:")
print(df["Gene_ID"].isna().sum())

print("\nDuplicate Gene_ID values:")
print(df["Gene_ID"].duplicated().sum())

print("\nUnique Gene_ID values:")
print(df["Gene_ID"].nunique())

# Remove any completely empty rows (just in case)
df = df.dropna(how="all")

# Remove rows without a Gene_ID
df = df[df["Gene_ID"].notna()]

# ============================
# Import into SQLite
# ============================

conn = sqlite3.connect(db_file)

df.to_sql(
    table_name,
    conn,
    if_exists="replace",
    index=False
)

conn.close()

print("\n========== IMPORT COMPLETE ==========")
print("Table:", table_name)
print("Rows imported:", len(df))
print("Columns:", list(df.columns))