import pandas as pd
import numpy as np
import os

print(
    "\n========== STEP 28: CONSERVED RESPONSE "
    "ANALYSIS ==========\n"
)

# =========================================================
# INPUT AND OUTPUT FILES
# =========================================================

INPUT_FILE = "Final_Clustered_Matrix.csv"

BIOLOGY_FILE = (
    "Final_Biological_Cluster_Characterization.csv"
)

OUTPUT_DIR = "Conserved_Response_Analysis_Outputs"

CONSERVATION_OUTPUT = "Conserved_Response_Genes.csv"
DIVERGENT_OUTPUT = "Divergent_Response_Genes.csv"
CANDIDATE_OUTPUT = "Top_Conserved_Response_Candidates.csv"
CLUSTER_SUMMARY_OUTPUT = "Conserved_Response_By_Cluster.csv"

CLUSTER_COL = "Final_Cluster"

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)

# =========================================================
# IDENTICAL THRESHOLDS TO STEP 25b
# =========================================================

HIGH_THRESHOLD = 0.75
MODERATE_THRESHOLD = 0.25

# =========================================================
# SELECT COMPARABLE MEASUREMENTS
# =========================================================
#
# Colistin expression comparison:
#   Pseudomonas RNA-seq versus Klebsiella K56 colistin log2FC
#
# TnSeq comparison:
#   PAO1 persister Mean SI versus KPPR1 log2FC
#
# KPPR1 raw fold change is deliberately excluded.
#

PSEUDO_COLISTIN = "Colistin_RNAseq_Fold_Change"
KLEB_COLISTIN = "K56_vs_Colistin_Log2FC"

PSEUDO_TNSEQ = "PAO1_Persister_All_Mean_SI"
KLEB_TNSEQ = "KPPR1_TnSeq_Log2FC_Output_Input"

# Optional descriptive measurements retained in the output.
PSEUDO_TOBRAMYCIN = "Tobramycin_RNAseq_Fold_Change"
KLEB_COMBINATION = "K56_vs_Combination_Log2FC"

# =========================================================
# LOAD FILES
# =========================================================

df = pd.read_csv(
    INPUT_FILE,
    low_memory=False
)

biology = pd.read_csv(
    BIOLOGY_FILE,
    low_memory=False
)

print(f"Loaded matrix: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

print(f"\nLoaded biological summary: {BIOLOGY_FILE}")
print(f"Biological-summary rows: {biology.shape[0]}")

# =========================================================
# CHECK REQUIRED COLUMNS
# =========================================================

required_matrix_cols = [
    CLUSTER_COL,
    PSEUDO_COLISTIN,
    KLEB_COLISTIN,
    PSEUDO_TNSEQ,
    KLEB_TNSEQ
]

missing_matrix_cols = [
    column
    for column in required_matrix_cols
    if column not in df.columns
]

if missing_matrix_cols:
    raise ValueError(
        "Missing required columns from the final clustered matrix:\n"
        + "\n".join(
            f"  - {column}"
            for column in missing_matrix_cols
        )
    )

required_biology_cols = [
    CLUSTER_COL,
    "Biological_Fingerprint",
    "Outlier_Status"
]

missing_biology_cols = [
    column
    for column in required_biology_cols
    if column not in biology.columns
]

if missing_biology_cols:
    raise ValueError(
        "Missing required columns from the biological summary:\n"
        + "\n".join(
            f"  - {column}"
            for column in missing_biology_cols
        )
    )

# Keep only one row per cluster before merging.
biology_for_merge = (
    biology[required_biology_cols]
    .drop_duplicates(subset=[CLUSTER_COL])
)

df = df.merge(
    biology_for_merge,
    on=CLUSTER_COL,
    how="left",
    validate="many_to_one"
)

# Convert relevant measurements to numeric.
numeric_cols = [
    PSEUDO_COLISTIN,
    KLEB_COLISTIN,
    PSEUDO_TNSEQ,
    KLEB_TNSEQ,
    PSEUDO_TOBRAMYCIN,
    KLEB_COMBINATION
]

for column in numeric_cols:
    if column in df.columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce"
        )

# =========================================================
# SHARED CLASSIFICATION FUNCTIONS
# =========================================================

def classify_expression(value):
    """Use the exact expression thresholds from Step 25b."""

    if pd.isna(value):
        return "NOT AVAILABLE"

    if value >= HIGH_THRESHOLD:
        return "HIGH"

    if value >= MODERATE_THRESHOLD:
        return "MOD HIGH"

    if value <= -HIGH_THRESHOLD:
        return "LOW"

    if value <= -MODERATE_THRESHOLD:
        return "MOD LOW"

    return "NEUTRAL"


def classify_fitness(value):
    """Use the exact fitness thresholds from Step 25b."""

    if pd.isna(value):
        return "NOT AVAILABLE"

    if value >= HIGH_THRESHOLD:
        return "HIGH FITNESS"

    if value >= MODERATE_THRESHOLD:
        return "MOD HIGH FITNESS"

    if value <= -HIGH_THRESHOLD:
        return "LOW FITNESS"

    if value <= -MODERATE_THRESHOLD:
        return "MOD LOW FITNESS"

    return "NEUTRAL FITNESS"


def expression_direction(value):
    """
    Convert a continuous expression value into broad direction.

    HIGH and MOD HIGH are positive.
    LOW and MOD LOW are negative.
    """

    category = classify_expression(value)

    if category in {"HIGH", "MOD HIGH"}:
        return "Positive"

    if category in {"LOW", "MOD LOW"}:
        return "Negative"

    if category == "NEUTRAL":
        return "Neutral"

    return "Not measured"


def fitness_direction(value):
    """
    Convert a mutant-fitness value into broad direction.
    """

    category = classify_fitness(value)

    if category in {
        "HIGH FITNESS",
        "MOD HIGH FITNESS"
    }:
        return "Positive"

    if category in {
        "LOW FITNESS",
        "MOD LOW FITNESS"
    }:
        return "Negative"

    if category == "NEUTRAL FITNESS":
        return "Neutral"

    return "Not measured"


# =========================================================
# CONSERVATION FUNCTIONS
# =========================================================

def conserved_expression(a, b):
    """
    Compare expression direction between Pseudomonas and Klebsiella.
    """

    a_direction = expression_direction(a)
    b_direction = expression_direction(b)

    if (
        a_direction == "Not measured"
        or b_direction == "Not measured"
    ):
        return "Not measured"

    if (
        a_direction == "Positive"
        and b_direction == "Positive"
    ):
        return "Conserved increased"

    if (
        a_direction == "Negative"
        and b_direction == "Negative"
    ):
        return "Conserved decreased"

    if (
        a_direction == "Positive"
        and b_direction == "Negative"
    ):
        return (
            "Divergent: Pseudomonas increased / "
            "Klebsiella decreased"
        )

    if (
        a_direction == "Negative"
        and b_direction == "Positive"
    ):
        return (
            "Divergent: Pseudomonas decreased / "
            "Klebsiella increased"
        )

    return "Not conserved / neutral"


def conserved_fitness(a, b):
    """
    Compare mutant-fitness direction between PAO1 and KPPR1.
    """

    a_direction = fitness_direction(a)
    b_direction = fitness_direction(b)

    if (
        a_direction == "Not measured"
        or b_direction == "Not measured"
    ):
        return "Not measured"

    if (
        a_direction == "Positive"
        and b_direction == "Positive"
    ):
        return "Shared increased mutant fitness"

    if (
        a_direction == "Negative"
        and b_direction == "Negative"
    ):
        return "Shared decreased mutant fitness"

    if (
        a_direction == "Positive"
        and b_direction == "Negative"
    ):
        return (
            "Divergent mutant fitness: "
            "Pseudomonas increased / Klebsiella decreased"
        )

    if (
        a_direction == "Negative"
        and b_direction == "Positive"
    ):
        return (
            "Divergent mutant fitness: "
            "Pseudomonas decreased / Klebsiella increased"
        )

    return "Not conserved / neutral mutant fitness"


def get_gene_label(row):
    """Choose the best available human-readable gene label."""

    candidate_columns = [
        "Gene_Label",
        "Pseudomonas_Gene_Name",
        "Klebsiella_Gene_Name",
        "Pseudomonas_Original_Locus_Tag",
        "Klebsiella_Locus_Tag"
    ]

    for column in candidate_columns:
        if column not in row.index:
            continue

        value = row[column]

        if pd.notna(value) and str(value).strip():
            return str(value).strip()

    return "Unknown"


# =========================================================
# GENE-LEVEL CLASSES
# =========================================================

df["Pseudo_Colistin_Class"] = (
    df[PSEUDO_COLISTIN].apply(classify_expression)
)

df["Kleb_Colistin_Class"] = (
    df[KLEB_COLISTIN].apply(classify_expression)
)

df["Pseudo_TnSeq_Class"] = (
    df[PSEUDO_TNSEQ].apply(classify_fitness)
)

df["Kleb_TnSeq_Class"] = (
    df[KLEB_TNSEQ].apply(classify_fitness)
)

# =========================================================
# CONSERVATION STATUS
# =========================================================

df["Colistin_Conservation_Status"] = df.apply(
    lambda row:
        conserved_expression(
            row[PSEUDO_COLISTIN],
            row[KLEB_COLISTIN]
        ),
    axis=1
)

df["TnSeq_Conservation_Status"] = df.apply(
    lambda row:
        conserved_fitness(
            row[PSEUDO_TNSEQ],
            row[KLEB_TNSEQ]
        ),
    axis=1
)

# =========================================================
# FLAGS AND SCORES
# =========================================================

df["Colistin_Conserved_Flag"] = (
    df["Colistin_Conservation_Status"].isin(
        [
            "Conserved increased",
            "Conserved decreased"
        ]
    )
)

df["TnSeq_Conserved_Flag"] = (
    df["TnSeq_Conservation_Status"].isin(
        [
            "Shared increased mutant fitness",
            "Shared decreased mutant fitness"
        ]
    )
)

df["Colistin_Divergent_Flag"] = (
    df["Colistin_Conservation_Status"]
    .str.startswith("Divergent", na=False)
)

df["TnSeq_Divergent_Flag"] = (
    df["TnSeq_Conservation_Status"]
    .str.startswith("Divergent", na=False)
)

df["Divergent_Response_Flag"] = (
    df["Colistin_Divergent_Flag"]
    | df["TnSeq_Divergent_Flag"]
)

df["Conserved_Response_Score"] = (
    df["Colistin_Conserved_Flag"].astype(int)
    + df["TnSeq_Conserved_Flag"].astype(int)
)

# Count how many direct cross-species comparisons were available.
df["Comparable_Response_Count"] = (
    (
        df[PSEUDO_COLISTIN].notna()
        & df[KLEB_COLISTIN].notna()
    ).astype(int)
    +
    (
        df[PSEUDO_TNSEQ].notna()
        & df[KLEB_TNSEQ].notna()
    ).astype(int)
)

# Response strength uses only the four directly compared values.
strength_cols = [
    PSEUDO_COLISTIN,
    KLEB_COLISTIN,
    PSEUDO_TNSEQ,
    KLEB_TNSEQ
]

df["Absolute_Response_Strength"] = (
    df[strength_cols]
    .abs()
    .mean(axis=1, skipna=True)
)

# A gene is not labelled conserved unless at least one direct
# comparison was actually available.
df["Has_Comparable_Data"] = (
    df["Comparable_Response_Count"] > 0
)

# =========================================================
# REPRESENTATIVE LABEL
# =========================================================

df["Representative_Gene_Label"] = df.apply(
    get_gene_label,
    axis=1
)

# =========================================================
# CONSERVED GENES
# =========================================================

conserved_df = df[
    (df["Conserved_Response_Score"] > 0)
    & df["Has_Comparable_Data"]
].copy()

conserved_df = conserved_df.sort_values(
    [
        "Conserved_Response_Score",
        "Absolute_Response_Strength"
    ],
    ascending=[False, False]
)

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
    "Pseudo_Colistin_Class",
    KLEB_COLISTIN,
    "Kleb_Colistin_Class",
    "Colistin_Conservation_Status",

    PSEUDO_TNSEQ,
    "Pseudo_TnSeq_Class",
    KLEB_TNSEQ,
    "Kleb_TnSeq_Class",
    "TnSeq_Conservation_Status",

    "Comparable_Response_Count",
    "Conserved_Response_Score",
    "Absolute_Response_Strength",
    "Outlier_Status"
]

conserved_cols = [
    column
    for column in conserved_cols
    if column in conserved_df.columns
]

conserved_path = os.path.join(
    OUTPUT_DIR,
    CONSERVATION_OUTPUT
)

conserved_df[conserved_cols].to_csv(
    conserved_path,
    index=False
)

# =========================================================
# DIVERGENT GENES
# =========================================================

divergent_df = df[
    df["Divergent_Response_Flag"]
    & df["Has_Comparable_Data"]
].copy()

divergent_df = divergent_df.sort_values(
    "Absolute_Response_Strength",
    ascending=False
)

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
    "Pseudo_Colistin_Class",
    KLEB_COLISTIN,
    "Kleb_Colistin_Class",
    "Colistin_Conservation_Status",

    PSEUDO_TNSEQ,
    "Pseudo_TnSeq_Class",
    KLEB_TNSEQ,
    "Kleb_TnSeq_Class",
    "TnSeq_Conservation_Status",

    "Comparable_Response_Count",
    "Absolute_Response_Strength",
    "Outlier_Status"
]

# Retain optional contextual measurements, but do not use them
# to define conservation or divergence.
if PSEUDO_TOBRAMYCIN in divergent_df.columns:
    divergent_cols.append(PSEUDO_TOBRAMYCIN)

if KLEB_COMBINATION in divergent_df.columns:
    divergent_cols.append(KLEB_COMBINATION)

divergent_cols = [
    column
    for column in divergent_cols
    if column in divergent_df.columns
]

divergent_path = os.path.join(
    OUTPUT_DIR,
    DIVERGENT_OUTPUT
)

divergent_df[divergent_cols].to_csv(
    divergent_path,
    index=False
)

# =========================================================
# TOP CONSERVED CANDIDATES
# =========================================================

candidate_df = conserved_df.head(50).copy()

candidate_path = os.path.join(
    OUTPUT_DIR,
    CANDIDATE_OUTPUT
)

candidate_df[conserved_cols].to_csv(
    candidate_path,
    index=False
)

# =========================================================
# CLUSTER-LEVEL SUMMARY
# =========================================================

cluster_summary = (
    df.groupby(CLUSTER_COL)
    .agg(
        Cluster_Size=(
            CLUSTER_COL,
            "size"
        ),

        Genes_With_Comparable_Data=(
            "Has_Comparable_Data",
            "sum"
        ),

        Conserved_Response_Genes=(
            "Conserved_Response_Score",
            lambda values:
                (values > 0).sum()
        ),

        Divergent_Response_Genes=(
            "Divergent_Response_Flag",
            "sum"
        ),

        Colistin_Conserved_Genes=(
            "Colistin_Conserved_Flag",
            "sum"
        ),

        TnSeq_Conserved_Genes=(
            "TnSeq_Conserved_Flag",
            "sum"
        ),

        Colistin_Divergent_Genes=(
            "Colistin_Divergent_Flag",
            "sum"
        ),

        TnSeq_Divergent_Genes=(
            "TnSeq_Divergent_Flag",
            "sum"
        ),

        Mean_Conserved_Response_Score=(
            "Conserved_Response_Score",
            "mean"
        ),

        Mean_Absolute_Response_Strength=(
            "Absolute_Response_Strength",
            "mean"
        )
    )
    .reset_index()
)

cluster_summary = cluster_summary.merge(
    biology_for_merge,
    on=CLUSTER_COL,
    how="left",
    validate="one_to_one"
)

cluster_summary["Percent_Conserved_Response_Genes"] = (
    np.where(
        cluster_summary["Genes_With_Comparable_Data"] > 0,

        cluster_summary["Conserved_Response_Genes"]
        / cluster_summary["Genes_With_Comparable_Data"]
        * 100,

        np.nan
    )
)

cluster_summary["Percent_Divergent_Response_Genes"] = (
    np.where(
        cluster_summary["Genes_With_Comparable_Data"] > 0,

        cluster_summary["Divergent_Response_Genes"]
        / cluster_summary["Genes_With_Comparable_Data"]
        * 100,

        np.nan
    )
)

cluster_summary = cluster_summary.sort_values(
    [
        "Conserved_Response_Genes",
        "Percent_Conserved_Response_Genes"
    ],
    ascending=[False, False]
)

cluster_summary_path = os.path.join(
    OUTPUT_DIR,
    CLUSTER_SUMMARY_OUTPUT
)

cluster_summary.to_csv(
    cluster_summary_path,
    index=False
)

# =========================================================
# PRINT OUTPUTS
# =========================================================

print(f"\nSaved conserved genes: {conserved_path}")
print(f"Saved divergent genes: {divergent_path}")
print(f"Saved top candidates: {candidate_path}")
print(f"Saved cluster summary: {cluster_summary_path}")

print("\nConserved response overview:")
print(
    "Genes with at least one comparable cross-species response: "
    f"{int(df['Has_Comparable_Data'].sum())}"
)
print(
    f"Total genes with a conserved response: "
    f"{len(conserved_df)}"
)
print(
    f"Total genes with a divergent response: "
    f"{len(divergent_df)}"
)

print("\nConserved-response score distribution:")
print(
    df.loc[
        df["Has_Comparable_Data"],
        "Conserved_Response_Score"
    ]
    .value_counts()
    .sort_index()
    .to_string()
)

print("\nTop clusters by conserved response genes:")

cluster_display_cols = [
    CLUSTER_COL,
    "Cluster_Size",
    "Genes_With_Comparable_Data",
    "Conserved_Response_Genes",
    "Percent_Conserved_Response_Genes",
    "Divergent_Response_Genes",
    "Percent_Divergent_Response_Genes",
    "Biological_Fingerprint"
]

print(
    cluster_summary[
        cluster_display_cols
    ]
    .head(15)
    .to_string(index=False)
)

print("\nTop conserved response candidates:")

if not candidate_df.empty:
    candidate_display_cols = [
        "Representative_Gene_Label",
        CLUSTER_COL,
        "Colistin_Conservation_Status",
        "TnSeq_Conservation_Status",
        "Conserved_Response_Score",
        "Absolute_Response_Strength"
    ]

    print(
        candidate_df[
            candidate_display_cols
        ]
        .head(20)
        .to_string(index=False)
    )
else:
    print(
        "No conserved response candidates were found using "
        "the current thresholds."
    )

print(
    "\n========== STEP 28 COMPLETE =========="
)