import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"

conn = sqlite3.connect(DB_NAME)

tnseq = pd.read_sql_query("""
SELECT
    Gene_ID,
    Gene_Name
FROM Klebsiella_TnSeq
""", conn)

kppr1 = pd.read_sql_query("""
SELECT
    Locus_Tag,
    Symbol,
    Protein_Accession,
    Protein_Name
FROM KPPR1_Genome_Annotation
""", conn)

print("TnSeq examples:")
print(tnseq.head(10))

print("\nKPPR1 annotation examples:")
print(kppr1.head(10))

# Check direct match: TnSeq Gene_ID vs KPPR1 Locus_Tag
direct = set(tnseq["Gene_ID"].dropna()).intersection(
    set(kppr1["Locus_Tag"].dropna())
)

# Check TnSeq Gene_ID vs KPPR1 Symbol
symbol_match = set(tnseq["Gene_ID"].dropna()).intersection(
    set(kppr1["Symbol"].dropna())
)

# Check TnSeq Gene_Name vs KPPR1 Symbol
gene_name_match = set(tnseq["Gene_Name"].dropna()).intersection(
    set(kppr1["Symbol"].dropna())
)

print("\nDirect Gene_ID ↔ Locus_Tag matches:", len(direct))
print("Gene_ID ↔ Symbol matches:", len(symbol_match))
print("Gene_Name ↔ Symbol matches:", len(gene_name_match))

if gene_name_match:
    print("\nExample Gene_Name ↔ Symbol matches:")
    print(list(gene_name_match)[:20])

conn.close()