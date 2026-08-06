import pandas as pd
import numpy as np

print(
    "\n========== STEP 25b: FINAL BIOLOGICAL "
    "CLUSTER CHARACTERIZATION ==========\n"
)

# =========================================================
# INPUT AND OUTPUT FILES
# =========================================================

INPUT_FILE = "Final_Clustered_Matrix.csv"

FINAL_OUTPUT = "Final_Biological_Cluster_Characterization.csv"
DETAILED_OUTPUT = (
    "Final_Biological_Cluster_Characterization_DETAILED.csv"
)
GENE_CATALOGUE_OUTPUT = "Final_Cluster_Gene_Catalogue.csv"

CLUSTER_COL = "Final_Cluster"

# =========================================================
# SHARED CLASSIFICATION THRESHOLDS
# =========================================================
#
# These same values must also be used in Step 28.
#
# Positive:
#   >= 0.75       HIGH
#   0.25 to 0.75  MOD HIGH
#
# Neutral:
#   -0.25 to 0.25 NEUTRAL
#
# Negative:
#   -0.75 to -0.25 MOD LOW
#   <= -0.75        LOW
#

HIGH_THRESHOLD = 0.75
MODERATE_THRESHOLD = 0.25

# =========================================================
# LOAD DATA
# =========================================================

df = pd.read_csv(INPUT_FILE, low_memory=False)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

if CLUSTER_COL not in df.columns:
    raise ValueError(
        f"Required cluster column '{CLUSTER_COL}' was not found."
    )

# =========================================================
# BIOLOGICAL GROUPS
# =========================================================
#
# Important corrections:
#
# 1. KPPR1 uses log2FC only.
#    Fold change and log2FC represent the same underlying ratio
#    and must not be averaged together.
#
# 2. PAO1 persister uses Mean SI only.
#    Rep1, Rep2, Rep3 and their calculated mean should not all
#    be included because this double-counts the same experiment.
#
# 3. Columns grouped together should be on comparable scales.
#

GROUPS = {
    "Pseudomonas_Colistin": {
        "cols": [
            "Colistin_RNAseq_Fold_Change",
            "Colistin_Riboseq_Fold_Change"
        ],
        "type": "expression",
        "short": "Colistin"
    },

    "Pseudomonas_Tobramycin": {
        "cols": [
            "Tobramycin_RNAseq_Fold_Change",
            "Tobramycin_Riboseq_Fold_Change"
        ],
        "type": "expression",
        "short": "Tobramycin"
    },

    "Pseudomonas_Persister_TnSeq": {
        "cols": [
            "PAO1_Persister_All_Mean_SI"
        ],
        "type": "fitness",
        "short": "Persister"
    },

    "Klebsiella_K56_Response": {
        "cols": [
            "K56_vs_Colistin_Log2FC",
            "K56_vs_Combination_Log2FC",
            "Colistin_vs_Combination_Log2FC"
        ],
        "type": "expression",
        "short": "K56"
    },

    "Klebsiella_KPPR1_TnSeq": {
        "cols": [
            "KPPR1_TnSeq_Log2FC_Output_Input"
        ],
        "type": "fitness",
        "short": "KPPR1"
    }
}

ID_COLS = [
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "Gene_Label"
]

# =========================================================
# CHECK REQUIRED COLUMNS
# =========================================================

missing_columns = []

for group_name, information in GROUPS.items():
    for column in information["cols"]:
        if column not in df.columns:
            missing_columns.append(column)

if missing_columns:
    missing_columns = sorted(set(missing_columns))

    raise ValueError(
        "The following required omics columns are missing:\n"
        + "\n".join(f"  - {column}" for column in missing_columns)
    )

# Convert all omics measurements to numeric.
all_omics_cols = []

for information in GROUPS.values():
    all_omics_cols.extend(information["cols"])

all_omics_cols = list(dict.fromkeys(all_omics_cols))

for column in all_omics_cols:
    df[column] = pd.to_numeric(df[column], errors="coerce")

# =========================================================
# OPTIONAL SCALE WARNING
# =========================================================
#
# The classification system assumes that expression values are
# signed and centred around zero.
#
# If an expression column contains only non-negative values and
# has a neutral reference near 1, it may be an ordinary fold-change
# ratio rather than a signed/log-transformed value.
#

expression_cols = []

for information in GROUPS.values():
    if information["type"] == "expression":
        expression_cols.extend(information["cols"])

print("\nExpression-scale checks:")

for column in expression_cols:
    non_null = df[column].dropna()

    if non_null.empty:
        print(f"  WARNING: {column} contains no numeric values.")
        continue

    minimum = non_null.min()
    median = non_null.median()
    maximum = non_null.max()

    print(
        f"  {column}: "
        f"min={minimum:.4f}, "
        f"median={median:.4f}, "
        f"max={maximum:.4f}"
    )

    if minimum >= 0 and 0.75 <= median <= 1.25:
        print(
            "    WARNING: This column may contain ordinary fold-change "
            "ratios centred near 1 rather than signed values centred at 0."
        )

# =========================================================
# CLASSIFICATION FUNCTIONS
# =========================================================

def classify_expression(value):
    """Classify a signed, zero-centred expression response."""

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
    """Classify a signed mutant-fitness measurement."""

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


def sentence_status(category, module_name):
    """Convert a category into presentation-friendly wording."""

    wording = {
        "HIGH": f"high {module_name} response",
        "MOD HIGH": f"moderately increased {module_name} response",
        "LOW": f"low {module_name} response",
        "MOD LOW": f"moderately reduced {module_name} response",
        "NEUTRAL": f"neutral {module_name} response",

        "HIGH FITNESS": f"high {module_name} mutant fitness",
        "MOD HIGH FITNESS":
            f"moderately high {module_name} mutant fitness",
        "LOW FITNESS": f"low {module_name} mutant fitness",
        "MOD LOW FITNESS":
            f"moderately low {module_name} mutant fitness",
        "NEUTRAL FITNESS":
            f"neutral {module_name} mutant fitness",

        "NOT AVAILABLE": f"{module_name} data not available"
    }

    return wording.get(
        category,
        f"{str(category).lower()} {module_name}"
    )


# =========================================================
# CLUSTER-LEVEL GROUP STATISTICS
# =========================================================

def add_group_stats(summary, group_name, cols, group_type):
    """
    Calculate cluster-level statistics for one biological group.

    Each selected column is reshaped into long format. This is valid
    only because the columns within each group are intended to be on
    comparable scales.

    KPPR1 fold change has deliberately been excluded.
    """

    values = df[[CLUSTER_COL] + cols].copy()

    long_values = values.melt(
        id_vars=CLUSTER_COL,
        value_vars=cols,
        var_name="Feature",
        value_name="Value"
    )

    long_values["Value"] = pd.to_numeric(
        long_values["Value"],
        errors="coerce"
    )

    stats = long_values.groupby(CLUSTER_COL)["Value"].agg(
        mean="mean",
        median="median",
        std="std",
        min="min",
        max="max",
        measured_count="count"
    )

    summary[f"{group_name}_Mean"] = stats["mean"]
    summary[f"{group_name}_Median"] = stats["median"]
    summary[f"{group_name}_Std"] = stats["std"]
    summary[f"{group_name}_Min"] = stats["min"]
    summary[f"{group_name}_Max"] = stats["max"]
    summary[f"{group_name}_Measured_Count"] = (
        stats["measured_count"]
    )

    if group_type == "fitness":
        summary[f"{group_name}_Percent_High_Fitness"] = (
            long_values.groupby(CLUSTER_COL)["Value"].apply(
                lambda values:
                    (values >= HIGH_THRESHOLD).mean() * 100
            )
        )

        summary[f"{group_name}_Percent_Low_Fitness"] = (
            long_values.groupby(CLUSTER_COL)["Value"].apply(
                lambda values:
                    (values <= -HIGH_THRESHOLD).mean() * 100
            )
        )

        summary[f"{group_name}_Category"] = (
            summary[f"{group_name}_Mean"].apply(classify_fitness)
        )

    else:
        summary[f"{group_name}_Percent_Induced"] = (
            long_values.groupby(CLUSTER_COL)["Value"].apply(
                lambda values:
                    (values >= HIGH_THRESHOLD).mean() * 100
            )
        )

        summary[f"{group_name}_Percent_Repressed"] = (
            long_values.groupby(CLUSTER_COL)["Value"].apply(
                lambda values:
                    (values <= -HIGH_THRESHOLD).mean() * 100
            )
        )

        summary[f"{group_name}_Category"] = (
            summary[f"{group_name}_Mean"].apply(
                classify_expression
            )
        )

    return summary


# =========================================================
# FINGERPRINT AND INTERPRETATION
# =========================================================

def make_fingerprint(row):
    """
    Build a cluster-level biological fingerprint.

    Every part of this fingerprint describes the cluster mean,
    not every individual gene in that cluster.
    """

    parts = []

    for group_name, information in GROUPS.items():
        category = row[f"{group_name}_Category"]
        short_name = information["short"]

        parts.append(f"{category} {short_name}")

    if row["Outlier_Status"] == "OUTLIER":
        parts.append("OUTLIER")

    return " | ".join(parts)


def get_dominant_drivers(row):
    """Identify the biological groups with the largest mean magnitude."""

    driver_scores = []

    for group_name, information in GROUPS.items():
        mean_value = row[f"{group_name}_Mean"]

        if pd.notna(mean_value):
            driver_scores.append(
                (
                    information["short"],
                    abs(mean_value),
                    mean_value
                )
            )

    driver_scores = sorted(
        driver_scores,
        key=lambda item: item[1],
        reverse=True
    )

    if len(driver_scores) == 0:
        return pd.Series({
            "Dominant_Driver": "Not available",
            "Dominant_Driver_Mean_Value": np.nan,
            "Secondary_Driver": "Not available",
            "Secondary_Driver_Mean_Value": np.nan
        })

    dominant = driver_scores[0]

    if len(driver_scores) > 1:
        secondary = driver_scores[1]
    else:
        secondary = ("Not available", np.nan, np.nan)

    return pd.Series({
        "Dominant_Driver": dominant[0],
        "Dominant_Driver_Mean_Value": dominant[2],
        "Secondary_Driver": secondary[0],
        "Secondary_Driver_Mean_Value": secondary[2]
    })


def make_biological_story(row):
    """Create a full-text description of the cluster profile."""

    module_names = {
        "Pseudomonas_Colistin":
            "Pseudomonas colistin",
        "Pseudomonas_Tobramycin":
            "Pseudomonas tobramycin",
        "Pseudomonas_Persister_TnSeq":
            "Pseudomonas persister TnSeq",
        "Klebsiella_K56_Response":
            "Klebsiella K56 antibiotic",
        "Klebsiella_KPPR1_TnSeq":
            "Klebsiella KPPR1 TnSeq"
    }

    story_parts = []

    for group_name in GROUPS:
        category = row[f"{group_name}_Category"]

        story_parts.append(
            sentence_status(
                category,
                module_names[group_name]
            )
        )

    story = ", ".join(story_parts) + "."

    if row["Outlier_Status"] == "OUTLIER":
        story += (
            " This cluster was also identified as an "
            "outlier cluster."
        )

    return story


def make_interpretation(row):
    """Generate a cautious cluster-level biological interpretation."""

    colistin = row["Pseudomonas_Colistin_Category"]
    tobramycin = row["Pseudomonas_Tobramycin_Category"]
    persister = row["Pseudomonas_Persister_TnSeq_Category"]
    k56 = row["Klebsiella_K56_Response_Category"]
    kppr1 = row["Klebsiella_KPPR1_TnSeq_Category"]

    positive_expression = {"HIGH", "MOD HIGH"}
    negative_expression = {"LOW", "MOD LOW"}

    positive_fitness = {
        "HIGH FITNESS",
        "MOD HIGH FITNESS"
    }

    negative_fitness = {
        "LOW FITNESS",
        "MOD LOW FITNESS"
    }

    if (
        colistin in positive_expression
        and k56 in positive_expression
    ):
        return (
            "This cluster shows increased antibiotic-response "
            "signals in both Pseudomonas and Klebsiella."
        )

    if (
        colistin in negative_expression
        and k56 in negative_expression
    ):
        return (
            "This cluster shows reduced antibiotic-response "
            "signals in both Pseudomonas and Klebsiella."
        )

    if colistin in positive_expression:
        return (
            "This cluster is mainly associated with an increased "
            "Pseudomonas colistin response."
        )

    if tobramycin in positive_expression:
        return (
            "This cluster is mainly associated with an increased "
            "Pseudomonas tobramycin response."
        )

    if k56 in positive_expression:
        return (
            "This cluster is mainly associated with an increased "
            "Klebsiella K56 antibiotic response."
        )

    if persister in negative_fitness:
        return (
            "This cluster contains mutants with reduced relative "
            "abundance in the Pseudomonas persister experiment."
        )

    if kppr1 in negative_fitness:
        return (
            "This cluster contains mutants with reduced relative "
            "abundance in the Klebsiella KPPR1 experiment."
        )

    if persister in positive_fitness:
        return (
            "This cluster contains mutants with increased relative "
            "abundance in the Pseudomonas persister experiment."
        )

    if kppr1 in positive_fitness:
        return (
            "This cluster contains mutants with increased relative "
            "abundance in the Klebsiella KPPR1 experiment."
        )

    return (
        "This cluster shows a mostly neutral or weak "
        "multi-omics response profile."
    )


# =========================================================
# REPRESENTATIVE GENES
# =========================================================

def representative_genes(cluster_id, top_n=10):
    """Select genes with the largest mean absolute response."""

    cluster_df = df[
        df[CLUSTER_COL] == cluster_id
    ].copy()

    cluster_df["Representative_Score"] = (
        cluster_df[all_omics_cols]
        .abs()
        .mean(axis=1, skipna=True)
    )

    if "Gene_Label" in cluster_df.columns:
        gene_col = "Gene_Label"
    elif "Pseudomonas_Gene_Name" in cluster_df.columns:
        gene_col = "Pseudomonas_Gene_Name"
    else:
        return ""

    genes = (
        cluster_df
        .sort_values(
            "Representative_Score",
            ascending=False
        )[gene_col]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .head(top_n)
        .tolist()
    )

    return "; ".join(genes)


# =========================================================
# BUILD CLUSTER SUMMARY
# =========================================================

summary = (
    df.groupby(CLUSTER_COL)
    .size()
    .reset_index(name="Cluster_Size")
    .set_index(CLUSTER_COL)
)

for group_name, information in GROUPS.items():
    summary = add_group_stats(
        summary=summary,
        group_name=group_name,
        cols=information["cols"],
        group_type=information["type"]
    )

mean_cols = [
    f"{group_name}_Mean"
    for group_name in GROUPS
]

summary["Overall_Response_Strength"] = (
    summary[mean_cols]
    .abs()
    .mean(axis=1, skipna=True)
)

# Outliers are clusters whose overall response strength is
# at least two standard deviations above the mean.
overall_mean = summary["Overall_Response_Strength"].mean()
overall_std = summary["Overall_Response_Strength"].std()

outlier_cutoff = overall_mean + (2 * overall_std)

summary["Outlier_Status"] = np.where(
    summary["Overall_Response_Strength"] >= outlier_cutoff,
    "OUTLIER",
    "NOT OUTLIER"
)

summary = summary.reset_index()

driver_df = summary.apply(
    get_dominant_drivers,
    axis=1
)

summary = pd.concat(
    [summary, driver_df],
    axis=1
)

summary["Biological_Fingerprint"] = summary.apply(
    make_fingerprint,
    axis=1
)

summary["Detailed_Biological_Story"] = summary.apply(
    make_biological_story,
    axis=1
)

summary["Biological_Interpretation"] = summary.apply(
    make_interpretation,
    axis=1
)

summary["Representative_Genes_Top10"] = (
    summary[CLUSTER_COL].apply(
        lambda cluster_id:
            representative_genes(cluster_id, top_n=10)
    )
)

summary = summary.sort_values(
    "Overall_Response_Strength",
    ascending=False
)

# =========================================================
# SAVE DETAILED SUMMARY
# =========================================================

summary.to_csv(
    DETAILED_OUTPUT,
    index=False
)

# =========================================================
# SAVE CLEAN FINAL SUMMARY
# =========================================================

final_cols = [
    CLUSTER_COL,
    "Cluster_Size",
    "Biological_Fingerprint",
    "Detailed_Biological_Story",
    "Dominant_Driver",
    "Dominant_Driver_Mean_Value",
    "Secondary_Driver",
    "Secondary_Driver_Mean_Value",
    "Overall_Response_Strength",
    "Outlier_Status",
    "Representative_Genes_Top10",
    "Biological_Interpretation"
]

summary[final_cols].to_csv(
    FINAL_OUTPUT,
    index=False
)

# =========================================================
# SAVE GENE CATALOGUE
# =========================================================

gene_catalogue_cols = (
    ID_COLS
    + all_omics_cols
    + [
        "Measured_Omics_Count",
        CLUSTER_COL
    ]
)

gene_catalogue_cols = [
    column
    for column in gene_catalogue_cols
    if column in df.columns
]

gene_catalogue = df[
    gene_catalogue_cols
].copy()

sort_cols = [CLUSTER_COL]

if "Gene_Label" in gene_catalogue.columns:
    sort_cols.append("Gene_Label")

gene_catalogue = gene_catalogue.sort_values(
    sort_cols
)

gene_catalogue.to_csv(
    GENE_CATALOGUE_OUTPUT,
    index=False
)

# =========================================================
# PRINT RESULTS
# =========================================================

print(f"\nSaved final summary: {FINAL_OUTPUT}")
print(f"Saved detailed summary: {DETAILED_OUTPUT}")
print(f"Saved gene catalogue: {GENE_CATALOGUE_OUTPUT}")

print("\nClassification thresholds:")
print(f"  HIGH: >= {HIGH_THRESHOLD}")
print(
    f"  MOD HIGH: >= {MODERATE_THRESHOLD} "
    f"and < {HIGH_THRESHOLD}"
)
print(
    f"  NEUTRAL: > {-MODERATE_THRESHOLD} "
    f"and < {MODERATE_THRESHOLD}"
)
print(
    f"  MOD LOW: > {-HIGH_THRESHOLD} "
    f"and <= {-MODERATE_THRESHOLD}"
)
print(f"  LOW: <= {-HIGH_THRESHOLD}")

print("\nFinal biological cluster profiles:")

display_cols = [
    CLUSTER_COL,
    "Cluster_Size",
    "Biological_Fingerprint",
    "Dominant_Driver",
    "Secondary_Driver",
    "Outlier_Status"
]

print(
    summary[display_cols].to_string(index=False)
)

print(
    "\n========== STEP 25b COMPLETE =========="
)