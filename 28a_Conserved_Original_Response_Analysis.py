import pandas as pd
import numpy as np
import os

print("\n========== STEP 28: CONSERVED RESPONSE ANALYSIS ==========\n")

INPUT_FILE = "Final_Clustered_Matrix.csv"
BIOLOGY_FILE = "Final_Biological_Cluster_Characterization.csv"

OUTPUT_DIR = "Conserved_Response_Analysis_Outputs"

CONSERVATION_OUTPUT = "Conserved_Response_Genes.csv"
CLUSTER_SUMMARY_OUTPUT = "Conserved_Response_By_Cluster.csv"
DIVERGENT_OUTPUT = "Divergent_Response_Genes.csv"
CANDIDATE_OUTPUT = "Top_Conserved_Response_Candidates.csv"

CLUSTER_COL = "Final_Cluster"

os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(INPUT_FILE)
biology = pd.read_csv(BIOLOGY_FILE)

print(f"Loaded matrix: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# Add biological fingerprint
df = df.merge(
    biology[[CLUSTER_COL, "Biological_Fingerprint", "Outlier_Status"]],
    on=CLUSTER_COL,
    how="left"
)

# ---------------------------------------------------
# Define response columns
# ---------------------------------------------------

PSEUDO_COLISTIN = "Colistin_RNAseq_Fold_Change"
PSEUDO_TOBRAMYCIN = "Tobramycin_RNAseq_Fold_Change"

KLEB_COLISTIN = "K56_vs_Colistin_Log2FC"
KLEB_COMBINATION = "K56_vs_Combination_Log2FC"

PSEUDO_TNSEQ = "PAO1_Persister_All_Mean_SI"
KLEB_TNSEQ = "KPPR1_TnSeq_Log2FC_Output_Input"

required_cols = [
    PSEUDO_COLISTIN,
    PSEUDO_TOBRAMYCIN,
    KLEB_COLISTIN,
    KLEB_COMBINATION,
    PSEUDO_TNSEQ,
    KLEB_TNSEQ,
    CLUSTER_COL
]

missing = [c for c in required_cols if c not in df.columns]

if missing:
    raise ValueError(f"Missing required columns: {missing}")

# ---------------------------------------------------
# Helper functions
# ---------------------------------------------------

def classify_expression(value):
    if pd.isna(value):
        return "Not measured"
    elif value >= 1:
        return "Induced"
    elif value <= -1:
        return "Repressed"
    else:
        return "Neutral"


def classify_fitness(value):
    if pd.isna(value):
        return "Not measured"
    elif value >= 1:
        return "High fitness"
    elif value <= -1:
        return "Low fitness"
    else:
        return "Neutral fitness"


def conserved_direction(a, b):
    """
    Classifies whether two responses move in the same or opposite direction.
    """
    a_class = classify_expression(a)
    b_class = classify_expression(b)

    if a_class == "Induced" and b_class == "Induced":
        return "Conserved induced"
    elif a_class == "Repressed" and b_class == "Repressed":
        return "Conserved repressed"
    elif a_class == "Induced" and b_class == "Repressed":
        return "Divergent: Pseudomonas induced / Klebsiella repressed"
    elif a_class == "Repressed" and b_class == "Induced":
        return "Divergent: Pseudomonas repressed / Klebsiella induced"
    else:
        return "Not conserved / neutral"


def conserved_fitness(a, b):
    a_class = classify_fitness(a)
    b_class = classify_fitness(b)

    if a_class == "Low fitness" and b_class == "Low fitness":
        return "Shared low fitness"
    elif a_class == "High fitness" and b_class == "High fitness":
        return "Shared high fitness"
    elif a_class == "Low fitness" and b_class == "High fitness":
        return "Divergent fitness: Pseudomonas low / Klebsiella high"
    elif a_class == "High fitness" and b_class == "Low fitness":
        return "Divergent fitness: Pseudomonas high / Klebsiella low"
    else:
        return "Not conserved / neutral fitness"


def get_gene_label(row):
    for col in [
        "Gene_Label",
        "Pseudomonas_Gene_Name",
        "Klebsiella_Gene_Name",
        "Pseudomonas_Original_Locus_Tag",
        "Klebsiella_Locus_Tag"
    ]:
        if col in row and pd.notna(row[col]) and str(row[col]).strip() != "":
            return str(row[col])
    return "Unknown"


# ---------------------------------------------------
# Expression conservation
# ---------------------------------------------------

df["Pseudo_Colistin_Class"] = df[PSEUDO_COLISTIN].apply(classify_expression)
df["Kleb_Colistin_Class"] = df[KLEB_COLISTIN].apply(classify_expression)

df["Pseudo_Tobramycin_Class"] = df[PSEUDO_TOBRAMYCIN].apply(classify_expression)
df["Kleb_Combination_Class"] = df[KLEB_COMBINATION].apply(classify_expression)

df["Colistin_Conservation_Status"] = df.apply(
    lambda row: conserved_direction(row[PSEUDO_COLISTIN], row[KLEB_COLISTIN]),
    axis=1
)

df["PseudoTobra_vs_KlebCombination_Status"] = df.apply(
    lambda row: conserved_direction(row[PSEUDO_TOBRAMYCIN], row[KLEB_COMBINATION]),
    axis=1
)

# ---------------------------------------------------
# Fitness conservation
# ---------------------------------------------------

df["Pseudo_TnSeq_Class"] = df[PSEUDO_TNSEQ].apply(classify_fitness)
df["Kleb_TnSeq_Class"] = df[KLEB_TNSEQ].apply(classify_fitness)

df["TnSeq_Conservation_Status"] = df.apply(
    lambda row: conserved_fitness(row[PSEUDO_TNSEQ], row[KLEB_TNSEQ]),
    axis=1
)

# ---------------------------------------------------
# Overall conserved response score
# ---------------------------------------------------

df["Colistin_Conserved_Flag"] = df["Colistin_Conservation_Status"].isin([
    "Conserved induced",
    "Conserved repressed"
])

df["TnSeq_Conserved_Flag"] = df["TnSeq_Conservation_Status"].isin([
    "Shared low fitness",
    "Shared high fitness"
])

df["Divergent_Response_Flag"] = (
    df["Colistin_Conservation_Status"].str.contains("Divergent", na=False)
    | df["PseudoTobra_vs_KlebCombination_Status"].str.contains("Divergent", na=False)
    | df["TnSeq_Conservation_Status"].str.contains("Divergent", na=False)
)

df["Conserved_Response_Score"] = (
    df["Colistin_Conserved_Flag"].astype(int)
    + df["TnSeq_Conserved_Flag"].astype(int)
)

df["Absolute_Response_Strength"] = df[
    [
        PSEUDO_COLISTIN,
        PSEUDO_TOBRAMYCIN,
        KLEB_COLISTIN,
        KLEB_COMBINATION,
        PSEUDO_TNSEQ,
        KLEB_TNSEQ
    ]
].abs().mean(axis=1)

# ---------------------------------------------------
# Save conserved genes
# ---------------------------------------------------

conserved_df = df[
    (df["Conserved_Response_Score"] > 0)
].copy()

conserved_df["Representative_Gene_Label"] = conserved_df.apply(get_gene_label, axis=1)

conserved_cols = [
    "Representative_Gene_Label",
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "Gene_Label",
    CLUSTER_COL,
    "Biological_Fingerprint",
    PSEUDO_COLISTIN,
    KLEB_COLISTIN,
    "Colistin_Conservation_Status",
    PSEUDO_TNSEQ,
    KLEB_TNSEQ,
    "TnSeq_Conservation_Status",
    "Conserved_Response_Score",
    "Absolute_Response_Strength",
    "Outlier_Status"
]

conserved_cols = [c for c in conserved_cols if c in conserved_df.columns]

conserved_df = conserved_df.sort_values(
    ["Conserved_Response_Score", "Absolute_Response_Strength"],
    ascending=[False, False]
)

conserved_path = os.path.join(OUTPUT_DIR, CONSERVATION_OUTPUT)
conserved_df[conserved_cols].to_csv(conserved_path, index=False)

# ---------------------------------------------------
# Save divergent genes
# ---------------------------------------------------

divergent_df = df[df["Divergent_Response_Flag"]].copy()
divergent_df["Representative_Gene_Label"] = divergent_df.apply(get_gene_label, axis=1)

divergent_cols = [
    "Representative_Gene_Label",
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "Gene_Label",
    CLUSTER_COL,
    "Biological_Fingerprint",
    PSEUDO_COLISTIN,
    KLEB_COLISTIN,
    "Colistin_Conservation_Status",
    PSEUDO_TOBRAMYCIN,
    KLEB_COMBINATION,
    "PseudoTobra_vs_KlebCombination_Status",
    PSEUDO_TNSEQ,
    KLEB_TNSEQ,
    "TnSeq_Conservation_Status",
    "Absolute_Response_Strength",
    "Outlier_Status"
]

divergent_cols = [c for c in divergent_cols if c in divergent_df.columns]

divergent_df = divergent_df.sort_values(
    "Absolute_Response_Strength",
    ascending=False
)

divergent_path = os.path.join(OUTPUT_DIR, DIVERGENT_OUTPUT)
divergent_df[divergent_cols].to_csv(divergent_path, index=False)

# ---------------------------------------------------
# Top conserved candidates
# ---------------------------------------------------

candidate_df = conserved_df.head(50).copy()

candidate_path = os.path.join(OUTPUT_DIR, CANDIDATE_OUTPUT)
candidate_df[conserved_cols].to_csv(candidate_path, index=False)

# ---------------------------------------------------
# Cluster-level conserved response summary
# ---------------------------------------------------

cluster_summary = df.groupby(CLUSTER_COL).agg(
    Cluster_Size=("Gene_Label", "count"),
    Conserved_Response_Genes=("Conserved_Response_Score", lambda x: (x > 0).sum()),
    Divergent_Response_Genes=("Divergent_Response_Flag", "sum"),
    Mean_Conserved_Response_Score=("Conserved_Response_Score", "mean"),
    Mean_Absolute_Response_Strength=("Absolute_Response_Strength", "mean")
).reset_index()

cluster_summary = cluster_summary.merge(
    biology[[CLUSTER_COL, "Biological_Fingerprint", "Outlier_Status"]],
    on=CLUSTER_COL,
    how="left"
)

cluster_summary["Percent_Conserved_Response_Genes"] = (
    cluster_summary["Conserved_Response_Genes"] / cluster_summary["Cluster_Size"] * 100
)

cluster_summary["Percent_Divergent_Response_Genes"] = (
    cluster_summary["Divergent_Response_Genes"] / cluster_summary["Cluster_Size"] * 100
)

cluster_summary = cluster_summary.sort_values(
    ["Conserved_Response_Genes", "Percent_Conserved_Response_Genes"],
    ascending=False
)

cluster_summary_path = os.path.join(OUTPUT_DIR, CLUSTER_SUMMARY_OUTPUT)
cluster_summary.to_csv(cluster_summary_path, index=False)

# ---------------------------------------------------
# Print outputs
# ---------------------------------------------------

print(f"\nSaved conserved genes: {conserved_path}")
print(f"Saved divergent genes: {divergent_path}")
print(f"Saved top candidates: {candidate_path}")
print(f"Saved cluster summary: {cluster_summary_path}")

print("\nConserved response overview:")
print(f"Total genes with conserved response: {len(conserved_df)}")
print(f"Total genes with divergent response: {len(divergent_df)}")

print("\nTop clusters by conserved response genes:")
print(
    cluster_summary[
        [
            CLUSTER_COL,
            "Cluster_Size",
            "Conserved_Response_Genes",
            "Percent_Conserved_Response_Genes",
            "Divergent_Response_Genes",
            "Percent_Divergent_Response_Genes",
            "Biological_Fingerprint"
        ]
    ].head(15).to_string(index=False)
)

print("\nTop conserved response candidates:")
if len(candidate_df) > 0:
    print(
        candidate_df[
            [
                "Representative_Gene_Label",
                CLUSTER_COL,
                "Colistin_Conservation_Status",
                "TnSeq_Conservation_Status",
                "Conserved_Response_Score",
                "Absolute_Response_Strength"
            ]
        ].head(20).to_string(index=False)
    )
else:
    print("No conserved response candidates found using the current thresholds.")

print("\n========== STEP 28a COMPLETE ==========")