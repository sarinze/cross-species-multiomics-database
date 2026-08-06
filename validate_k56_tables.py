import pandas as pd
import sqlite3

conn = sqlite3.connect("pseudomonas_new.db")

# Load tables
raw = pd.read_sql(
    "SELECT Gene_ID FROM KPHS_K56_COLISTIN_COMBINATION_RAW_COUNTS",
    conn
)

deg = pd.read_sql(
    "SELECT Gene FROM KPHS_K56_COLISTIN_COMBINATION_DEG",
    conn
)

conn.close()

# Convert to sets
raw_genes = set(raw["Gene_ID"])
deg_genes = set(deg["Gene"])

# Compare
shared = raw_genes & deg_genes
raw_only = raw_genes - deg_genes
deg_only = deg_genes - raw_genes

print("Raw genes:", len(raw_genes))
print("DEG genes:", len(deg_genes))
print("Shared genes:", len(shared))
print("Raw only:", len(raw_only))
print("DEG only:", len(deg_only))

if len(raw_only) > 0:
    print("\nExamples in RAW only:")
    print(list(raw_only)[:10])

if len(deg_only) > 0:
    print("\nExamples in DEG only:")
    print(list(deg_only)[:10])