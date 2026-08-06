import sqlite3
import pandas as pd

DB_NAME = "pseudomonas_new.db"

conn = sqlite3.connect(DB_NAME)

tnseq = pd.read_sql_query("""
SELECT
    Gene_ID,
    Gene_Name
FROM Klebsiella_TnSeq
WHERE Gene_ID IS NOT NULL
ORDER BY Gene_ID
""", conn)

kppr1 = pd.read_sql_query("""
SELECT
    Locus_Tag,
    Start,
    End,
    Orientation,
    Symbol,
    Protein_Name,
    Protein_Accession
FROM KPPR1_Genome_Annotation
WHERE Protein_Accession IS NOT NULL
ORDER BY Start
""", conn)

# Extract numeric part
tnseq["Old_Number"] = tnseq["Gene_ID"].str.extract(r"VK055_(\d+)").astype(float)
kppr1["RS_Number"] = kppr1["Locus_Tag"].str.extract(r"VK055_RS(\d+)").astype(float)

# Sort both by genome/order number
tnseq_sorted = tnseq.sort_values("Old_Number").reset_index(drop=True)
kppr1_sorted = kppr1.sort_values("Start").reset_index(drop=True)

# Create tentative row-order mapping
n = min(len(tnseq_sorted), len(kppr1_sorted))

mapping_test = pd.concat(
    [
        tnseq_sorted.iloc[:n].reset_index(drop=True),
        kppr1_sorted.iloc[:n].reset_index(drop=True)
    ],
    axis=1
)

# Check whether gene names match annotation symbols or protein names
mapping_test["GeneName_Equals_Symbol"] = (
    mapping_test["Gene_Name"].astype(str).str.lower()
    ==
    mapping_test["Symbol"].astype(str).str.lower()
)

mapping_test["GeneName_In_ProteinName"] = mapping_test.apply(
    lambda row: str(row["Gene_Name"]).lower() in str(row["Protein_Name"]).lower(),
    axis=1
)

print("TnSeq genes:", len(tnseq_sorted))
print("KPPR1 protein-coding annotation rows:", len(kppr1_sorted))
print("Tentative mapped rows:", len(mapping_test))

print("\nFirst 30 tentative mappings:")
print(mapping_test[
    [
        "Gene_ID",
        "Gene_Name",
        "Locus_Tag",
        "Symbol",
        "Protein_Name",
        "GeneName_Equals_Symbol",
        "GeneName_In_ProteinName"
    ]
].head(30))

print("\nValidation counts:")
print("GeneName equals Symbol:", mapping_test["GeneName_Equals_Symbol"].sum())
print("GeneName appears in ProteinName:", mapping_test["GeneName_In_ProteinName"].sum())

conn.close()