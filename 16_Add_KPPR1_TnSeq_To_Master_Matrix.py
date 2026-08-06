import sqlite3
import pandas as pd
import re

DB = "pseudomonas_new.db"

MASTER_INPUT = "Master_Matrix_With_Klebsiella_DEG.csv"
MASTER_OUTPUT = "Master_Multiomics_Matrix_With_KPPR1_TnSeq.csv"
MAPPING_CHECK_OUTPUT = "KPPR1_TnSeq_Mapping_Check.csv"

TNSEQ_TABLE = "Klebsiella_TnSeq"
LINKER_TABLE = "KPPR1_MGH78578_ANNOTATED_LINKER"
SQLITE_OUTPUT_TABLE = "Master_Multiomics_Matrix_With_KPPR1_TnSeq"

MASTER_KLEB_ID = "Klebsiella_Locus_Tag"
TNSEQ_KPPR1_ID = "Gene_ID"

LINKER_KPPR1_LOCUS = "KPPR1_Locus_Tag"
LINKER_MGH_LOCUS = "MGH78578_Locus_Tag"


def convert_vk055_to_rs(gene_id):
    """
    Converts old KPPR1 IDs:
    VK055_0002 -> VK055_RS00010
    VK055_0003 -> VK055_RS00015
    """
    if pd.isna(gene_id):
        return None

    gene_id = str(gene_id).strip()

    match = re.match(r"VK055_(\d+)$", gene_id)
    if not match:
        return gene_id

    number = int(match.group(1))
    rs_number = number * 5

    return f"VK055_RS{rs_number:05d}"


conn = sqlite3.connect(DB)

master = pd.read_csv(MASTER_INPUT, dtype=str, low_memory=False)
tnseq = pd.read_sql(f"SELECT * FROM {TNSEQ_TABLE}", conn)
linker = pd.read_sql(f"SELECT * FROM {LINKER_TABLE}", conn)

print("========== INPUT CHECK ==========")
print("Master rows:", len(master))
print("TnSeq rows:", len(tnseq))
print("Linker rows:", len(linker))

# Clean IDs
master[MASTER_KLEB_ID] = master[MASTER_KLEB_ID].astype(str).str.strip()
tnseq[TNSEQ_KPPR1_ID] = tnseq[TNSEQ_KPPR1_ID].astype(str).str.strip()
linker[LINKER_KPPR1_LOCUS] = linker[LINKER_KPPR1_LOCUS].astype(str).str.strip()
linker[LINKER_MGH_LOCUS] = linker[LINKER_MGH_LOCUS].astype(str).str.strip()

# Convert VK055_ IDs to VK055_RS IDs
tnseq["KPPR1_Converted_RS_Locus_Tag"] = tnseq[TNSEQ_KPPR1_ID].apply(convert_vk055_to_rs)

print("\n========== ID CONVERSION CHECK ==========")
print(
    tnseq[
        [TNSEQ_KPPR1_ID, "KPPR1_Converted_RS_Locus_Tag"]
    ].head(20)
)

# Rename TnSeq columns
tnseq_renamed = tnseq.rename(
    columns={
        "Gene_Name": "KPPR1_TnSeq_Gene_Name",
        "Gene_Length": "KPPR1_TnSeq_Gene_Length",
        "Total_Insertions": "KPPR1_TnSeq_Total_Insertions",
        "Unique_Insertions": "KPPR1_TnSeq_Unique_Insertions",
        "Mean_Input": "KPPR1_TnSeq_Mean_Input",
        "Mean_Output": "KPPR1_TnSeq_Mean_Output",
        "Fold_Change": "KPPR1_TnSeq_Fold_Change",
        "Log2_Fold_Change_Output_Input": "KPPR1_TnSeq_Log2FC_Output_Input",
        "P_Value": "KPPR1_TnSeq_P_Value",
        "Adjusted_P_Value": "KPPR1_TnSeq_Adjusted_P_Value",
        "Primary_KEGG_Annotation": "KPPR1_TnSeq_Primary_KEGG_Annotation",
        "Secondary_KEGG_Annotation": "KPPR1_TnSeq_Secondary_KEGG_Annotation",
    }
)

check_cols = [
    "KPPR1_TnSeq_Log2FC_Output_Input",
    "KPPR1_TnSeq_Adjusted_P_Value",
]

print("\n========== TNSEQ VALUE CHECK BEFORE MAPPING ==========")
print(tnseq_renamed[check_cols].notna().sum())

# Prepare linker
linker_keep_cols = [
    "KPPR1_Locus_Tag",
    "KPPR1_Protein_Accession",
    "KPPR1_Protein_Name",
    "KPPR1_Symbol",
    "MGH78578_Locus_Tag",
    "MGH78578_Protein_Accession",
    "MGH78578_Protein_Name",
    "MGH78578_Symbol",
    "pident_KPPR1_to_MGH78578",
    "bitscore_KPPR1_to_MGH78578",
]

linker_keep = linker[linker_keep_cols].drop_duplicates()

# Map TnSeq to linker
mapped_tnseq = tnseq_renamed.merge(
    linker_keep,
    left_on="KPPR1_Converted_RS_Locus_Tag",
    right_on="KPPR1_Locus_Tag",
    how="left",
)

print("\n========== KPPR1 TNSEQ MAPPING QC ==========")
print("TnSeq rows:", len(mapped_tnseq))
print("Mapped to MGH78578 locus:", mapped_tnseq["MGH78578_Locus_Tag"].notna().sum())
print("Unmapped:", mapped_tnseq["MGH78578_Locus_Tag"].isna().sum())

print("\nTnSeq values AFTER mapping:")
print(mapped_tnseq[check_cols].notna().sum())

print("\nMapped preview:")
print(
    mapped_tnseq[
        [
            "Gene_ID",
            "KPPR1_Converted_RS_Locus_Tag",
            "KPPR1_Locus_Tag",
            "MGH78578_Locus_Tag",
            "KPPR1_TnSeq_Log2FC_Output_Input",
            "KPPR1_TnSeq_Adjusted_P_Value",
        ]
    ].head(20)
)

# Check overlap with master
master_ids = set(master[MASTER_KLEB_ID].dropna().astype(str))
mapped_ids = set(mapped_tnseq["MGH78578_Locus_Tag"].dropna().astype(str))

overlap = master_ids & mapped_ids

print("\n========== MASTER OVERLAP CHECK ==========")
print("Master Klebsiella IDs:", len(master_ids))
print("Mapped MGH78578 IDs:", len(mapped_ids))
print("Overlap:", len(overlap))

# Save mapping check
mapped_tnseq.to_csv(MAPPING_CHECK_OUTPUT, index=False)

# Drop duplicate MGH mappings before merging
mapped_tnseq_for_master = mapped_tnseq.dropna(
    subset=["MGH78578_Locus_Tag"]
).drop_duplicates(
    subset=["MGH78578_Locus_Tag"]
)

# Merge into master
final = master.merge(
    mapped_tnseq_for_master,
    left_on=MASTER_KLEB_ID,
    right_on="MGH78578_Locus_Tag",
    how="left",
    suffixes=("", "_KPPR1Mapped"),
)

print("\n========== FINAL MERGE QC ==========")
print("Final rows:", len(final))
print("Final columns:", len(final.columns))

print("\nNon-null KPPR1 TnSeq values:")
print(final[check_cols].notna().sum())

print("\nMissing KPPR1 TnSeq values:")
print(final[check_cols].isna().sum())

# Save
final.to_csv(MASTER_OUTPUT, index=False)
final.to_sql(SQLITE_OUTPUT_TABLE, conn, if_exists="replace", index=False)

conn.close()

print("\n========== COMPLETE ==========")
print("Output CSV:", MASTER_OUTPUT)
print("SQLite table:", SQLITE_OUTPUT_TABLE)
print("Mapping check:", MAPPING_CHECK_OUTPUT)