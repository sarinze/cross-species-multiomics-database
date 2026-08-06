import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist


# ============================================================
# STEP 20B: HIERARCHICAL CLUSTERING USING ALL MASTER-MATRIX ROWS
# ============================================================

INPUT_FILE = "Master_Multiomics_Matrix_Normalized.csv"
NORMALIZED_COLUMNS_FILE = "Normalized_Columns.txt"

N_CLUSTERS = 20

OUTPUT_ASSIGNMENTS = "All_Genes_20_Cluster_Assignments.csv"
OUTPUT_CLUSTER_SUMMARY = "All_Genes_20_Cluster_Summary.csv"
OUTPUT_INFORMATION_SUMMARY = "All_Genes_Information_Content_Summary.csv"
OUTPUT_CLUSTER_INFORMATION = "All_Genes_Cluster_Information_Distribution.csv"
OUTPUT_PROFILE_MATRIX = "All_Genes_20_Cluster_Profile_Matrix.csv"
OUTPUT_DENDROGRAM = "All_Genes_20_Cluster_Dendrogram.png"
OUTPUT_LINKAGE = "All_Genes_Hierarchical_Linkage_Matrix.npy"


print("\n========== STEP 20B: CLUSTER ALL MASTER-MATRIX GENES ==========\n")


# ------------------------------------------------------------
# 1. Check that the required files exist
# ------------------------------------------------------------

for required_file in [INPUT_FILE, NORMALIZED_COLUMNS_FILE]:
    if not os.path.exists(required_file):
        raise FileNotFoundError(
            f"Required file not found: {required_file}\n"
            "Make sure the script is being run from the project folder."
        )


# ------------------------------------------------------------
# 2. Load the normalized master matrix
# ------------------------------------------------------------

df = pd.read_csv(INPUT_FILE, low_memory=False)

print(f"Loaded normalized matrix: {INPUT_FILE}")
print(f"Rows in complete master matrix: {df.shape[0]}")
print(f"Columns in complete master matrix: {df.shape[1]}")


# ------------------------------------------------------------
# 3. Read the normalized omics-column list
# ------------------------------------------------------------

with open(NORMALIZED_COLUMNS_FILE, "r", encoding="utf-8") as file:
    omics_cols = [line.strip() for line in file if line.strip()]

missing_omics_cols = [col for col in omics_cols if col not in df.columns]

if missing_omics_cols:
    raise ValueError(
        "The following normalized columns listed in "
        f"{NORMALIZED_COLUMNS_FILE} were not found in {INPUT_FILE}:\n"
        + "\n".join(missing_omics_cols)
    )

if len(omics_cols) < 2:
    raise ValueError("At least two omics columns are required for clustering.")

print(f"Normalized omics columns used: {len(omics_cols)}")

for col in omics_cols:
    print(f"  - {col}")


# ------------------------------------------------------------
# 4. Convert the omics block to numeric values
# ------------------------------------------------------------

omics_numeric = df[omics_cols].apply(
    pd.to_numeric,
    errors="coerce"
)

df[omics_cols] = omics_numeric


# ------------------------------------------------------------
# 5. Measure the information available for every row
# ------------------------------------------------------------

df["Measured_Omics_Count"] = omics_numeric.notna().sum(axis=1)
df["Missing_Omics_Count"] = len(omics_cols) - df["Measured_Omics_Count"]
df["Measured_Omics_Percent"] = (
    df["Measured_Omics_Count"] / len(omics_cols) * 100
).round(2)


def assign_information_level(measured_count):
    """
    Categorise each row according to how many of the 16 omics
    measurements are present.
    """
    if measured_count <= 2:
        return "Low (0-2 measured)"
    elif measured_count <= 5:
        return "Moderate (3-5 measured)"
    elif measured_count <= 10:
        return "High (6-10 measured)"
    else:
        return "Very high (11-16 measured)"


df["Information_Level"] = df["Measured_Omics_Count"].apply(
    assign_information_level
)


# ------------------------------------------------------------
# 6. Create readable gene/ortholog labels
# ------------------------------------------------------------

label_source_columns = [
    "Ortholog_Pair_ID",
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "KPPR1_Locus_Tag",
    "KPPR1_TnSeq_Gene_Name"
]


def safe_text(row, column):
    if column not in row.index:
        return ""

    value = row[column]

    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.lower() in ["nan", "none"]:
        return ""

    return value


def build_gene_label(row):
    values = []

    for column in label_source_columns:
        value = safe_text(row, column)

        if value and value not in values:
            values.append(value)

    if values:
        return " | ".join(values)

    return f"Master_Row_{row.name + 1}"


df["Gene_Label"] = df.apply(build_gene_label, axis=1)


# ------------------------------------------------------------
# 7. Create the all-gene clustering matrix
# ------------------------------------------------------------

# The normalized master matrix itself is not changed.
#
# For clustering only, NaN is filled with zero.
# In a z-score-scaled column, zero represents the column mean.
#
# This permits every master-matrix row to enter the exploratory
# clustering. Rows with little or no measured data are retained and
# labelled, rather than excluded.

X = omics_numeric.fillna(0.0)

if X.isna().any().any():
    raise ValueError(
        "Missing values remain in the clustering matrix after filling."
    )

if np.isinf(X.to_numpy()).any():
    raise ValueError(
        "Infinite values were detected in the clustering matrix."
    )

print("\nInformation-content distribution before clustering:")

information_order = [
    "Low (0-2 measured)",
    "Moderate (3-5 measured)",
    "High (6-10 measured)",
    "Very high (11-16 measured)"
]

information_counts = (
    df["Information_Level"]
    .value_counts()
    .reindex(information_order, fill_value=0)
)

for information_level, count in information_counts.items():
    percent = count / len(df) * 100

    print(
        f"  {information_level}: "
        f"{count} rows ({percent:.2f}%)"
    )

zero_measured_count = int((df["Measured_Omics_Count"] == 0).sum())
one_or_two_count = int(
    df["Measured_Omics_Count"].between(1, 2).sum()
)

print(f"\nRows with 0 measured omics values: {zero_measured_count}")
print(f"Rows with only 1-2 measured omics values: {one_or_two_count}")
print(f"Rows entering clustering: {X.shape[0]}")
print(f"Variables entering clustering: {X.shape[1]}")


# ------------------------------------------------------------
# 8. Perform hierarchical clustering
# ------------------------------------------------------------

if X.shape[0] < 2:
    raise ValueError("At least two rows are required for clustering.")

print("\nCalculating Euclidean distances for all genes...")

distance_matrix = pdist(
    X.to_numpy(dtype=np.float64),
    metric="euclidean"
)

print("Constructing Ward hierarchical linkage...")

linkage_matrix = linkage(
    distance_matrix,
    method="ward"
)

np.save(
    OUTPUT_LINKAGE,
    linkage_matrix
)

print(f"Saved linkage matrix: {OUTPUT_LINKAGE}")


# ------------------------------------------------------------
# 9. Assign 20 clusters
# ------------------------------------------------------------

cluster_labels = fcluster(
    linkage_matrix,
    t=N_CLUSTERS,
    criterion="maxclust"
)

df["All_Genes_Cluster_20"] = cluster_labels

actual_clusters = df["All_Genes_Cluster_20"].nunique()

print(f"\nClusters requested: {N_CLUSTERS}")
print(f"Clusters produced: {actual_clusters}")


# ------------------------------------------------------------
# 10. Save the complete assignment table
# ------------------------------------------------------------

preferred_metadata_columns = [
    "Ortholog_Pair_ID",
    "Pseudomonas_Strain",
    "Klebsiella_Strain",
    "Pseudomonas_Locus_Tag",
    "Pseudomonas_Original_Locus_Tag",
    "Pseudomonas_Gene_Name",
    "Pseudomonas_Symbol",
    "Klebsiella_Locus_Tag",
    "Klebsiella_Gene_Name",
    "Klebsiella_Symbol",
    "Gene_ID",
    "KPPR1_TnSeq_Gene_Name",
    "KPPR1_Converted_RS_Locus_Tag",
    "KPPR1_Locus_Tag",
    "KPPR1_Protein_Name",
    "KPPR1_Symbol",
    "MGH78578_Locus_Tag",
    "MGH78578_Protein_Name",
    "MGH78578_Symbol",
    "Gene_Label",
    "Measured_Omics_Count",
    "Missing_Omics_Count",
    "Measured_Omics_Percent",
    "Information_Level",
    "All_Genes_Cluster_20"
]

available_metadata_columns = [
    col for col in preferred_metadata_columns
    if col in df.columns
]

assignment_columns = (
    available_metadata_columns
    + [
        col for col in omics_cols
        if col not in available_metadata_columns
    ]
)

assignment_output = df[assignment_columns].copy()

assignment_output.to_csv(
    OUTPUT_ASSIGNMENTS,
    index=False
)

print(f"Saved complete assignments: {OUTPUT_ASSIGNMENTS}")


# ------------------------------------------------------------
# 11. Create the overall information-content summary
# ------------------------------------------------------------

information_summary = (
    df.groupby(
        "Information_Level",
        observed=False
    )
    .agg(
        Gene_Count=("Gene_Label", "size"),
        Mean_Measured_Omics=("Measured_Omics_Count", "mean"),
        Median_Measured_Omics=("Measured_Omics_Count", "median"),
        Minimum_Measured_Omics=("Measured_Omics_Count", "min"),
        Maximum_Measured_Omics=("Measured_Omics_Count", "max")
    )
    .reset_index()
)

information_summary["Percent_of_All_Rows"] = (
    information_summary["Gene_Count"] / len(df) * 100
).round(2)

information_summary["Information_Level"] = pd.Categorical(
    information_summary["Information_Level"],
    categories=information_order,
    ordered=True
)

information_summary = information_summary.sort_values(
    "Information_Level"
)

information_summary.to_csv(
    OUTPUT_INFORMATION_SUMMARY,
    index=False
)

print(
    "Saved information-content summary: "
    f"{OUTPUT_INFORMATION_SUMMARY}"
)


# ------------------------------------------------------------
# 12. Create cluster-level summaries
# ------------------------------------------------------------

cluster_summary = (
    df.groupby("All_Genes_Cluster_20")
    .agg(
        Gene_Count=("Gene_Label", "size"),
        Mean_Measured_Omics_Count=("Measured_Omics_Count", "mean"),
        Median_Measured_Omics_Count=("Measured_Omics_Count", "median"),
        Minimum_Measured_Omics_Count=("Measured_Omics_Count", "min"),
        Maximum_Measured_Omics_Count=("Measured_Omics_Count", "max"),
        Rows_With_Zero_Measurements=(
            "Measured_Omics_Count",
            lambda series: int((series == 0).sum())
        ),
        Rows_With_One_or_Two_Measurements=(
            "Measured_Omics_Count",
            lambda series: int(series.between(1, 2).sum())
        ),
        Rows_With_At_Least_Six_Measurements=(
            "Measured_Omics_Count",
            lambda series: int((series >= 6).sum())
        )
    )
    .reset_index()
)

cluster_summary["Percent_of_All_Rows"] = (
    cluster_summary["Gene_Count"] / len(df) * 100
).round(2)

cluster_summary["Percent_Zero_Measurements"] = (
    cluster_summary["Rows_With_Zero_Measurements"]
    / cluster_summary["Gene_Count"]
    * 100
).round(2)

cluster_summary["Percent_With_At_Least_Six_Measurements"] = (
    cluster_summary["Rows_With_At_Least_Six_Measurements"]
    / cluster_summary["Gene_Count"]
    * 100
).round(2)


# ------------------------------------------------------------
# 13. Calculate cluster means for the omics variables
# ------------------------------------------------------------

# These means use the clustering representation, where missing values
# are represented as zero. They therefore describe the profiles used
# by the clustering algorithm.

filled_profile_data = X.copy()
filled_profile_data["All_Genes_Cluster_20"] = cluster_labels

cluster_profile_filled = (
    filled_profile_data
    .groupby("All_Genes_Cluster_20")[omics_cols]
    .mean()
    .reset_index()
)

cluster_profile_filled.columns = [
    "All_Genes_Cluster_20"
    if col == "All_Genes_Cluster_20"
    else f"Filled_Mean_{col}"
    for col in cluster_profile_filled.columns
]


# These means use only genuinely measured values and ignore missing
# values. They are more appropriate for biological interpretation.

measured_profile_data = omics_numeric.copy()
measured_profile_data["All_Genes_Cluster_20"] = cluster_labels

cluster_profile_measured = (
    measured_profile_data
    .groupby("All_Genes_Cluster_20")[omics_cols]
    .mean()
    .reset_index()
)

cluster_profile_measured.columns = [
    "All_Genes_Cluster_20"
    if col == "All_Genes_Cluster_20"
    else f"Measured_Only_Mean_{col}"
    for col in cluster_profile_measured.columns
]


cluster_summary = cluster_summary.merge(
    cluster_profile_measured,
    on="All_Genes_Cluster_20",
    how="left"
)

cluster_summary.to_csv(
    OUTPUT_CLUSTER_SUMMARY,
    index=False
)

print(f"Saved cluster summary: {OUTPUT_CLUSTER_SUMMARY}")


profile_matrix = cluster_profile_filled.merge(
    cluster_profile_measured,
    on="All_Genes_Cluster_20",
    how="left"
)

profile_matrix.to_csv(
    OUTPUT_PROFILE_MATRIX,
    index=False
)

print(f"Saved cluster profile matrix: {OUTPUT_PROFILE_MATRIX}")


# ------------------------------------------------------------
# 14. Show information levels inside each cluster
# ------------------------------------------------------------

cluster_information_counts = pd.crosstab(
    df["All_Genes_Cluster_20"],
    df["Information_Level"]
)

for information_level in information_order:
    if information_level not in cluster_information_counts.columns:
        cluster_information_counts[information_level] = 0

cluster_information_counts = cluster_information_counts[
    information_order
]

cluster_information_percent = (
    cluster_information_counts
    .div(cluster_information_counts.sum(axis=1), axis=0)
    .mul(100)
    .round(2)
)

cluster_information_counts.columns = [
    f"Count_{col}" for col in cluster_information_counts.columns
]

cluster_information_percent.columns = [
    f"Percent_{col}" for col in cluster_information_percent.columns
]

cluster_information_output = pd.concat(
    [
        cluster_information_counts,
        cluster_information_percent
    ],
    axis=1
).reset_index()

cluster_information_output.to_csv(
    OUTPUT_CLUSTER_INFORMATION,
    index=False
)

print(
    "Saved cluster information distribution: "
    f"{OUTPUT_CLUSTER_INFORMATION}"
)


# ------------------------------------------------------------
# 15. Generate an all-gene dendrogram
# ------------------------------------------------------------

print("\nGenerating dendrogram...")

plt.figure(figsize=(16, 9))

dendrogram(
    linkage_matrix,
    no_labels=True,
    color_threshold=None
)

plt.title(
    "Hierarchical Clustering of All Master Multi-Omics Rows"
)
plt.xlabel("Genes / Ortholog Rows")
plt.ylabel("Ward Linkage Distance")

plt.tight_layout()

plt.savefig(
    OUTPUT_DENDROGRAM,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Saved dendrogram: {OUTPUT_DENDROGRAM}")


# ------------------------------------------------------------
# 16. Print the final cluster sizes and information content
# ------------------------------------------------------------

print("\nCluster sizes and data coverage:")

display_columns = [
    "All_Genes_Cluster_20",
    "Gene_Count",
    "Percent_of_All_Rows",
    "Mean_Measured_Omics_Count",
    "Rows_With_Zero_Measurements",
    "Rows_With_One_or_Two_Measurements",
    "Rows_With_At_Least_Six_Measurements",
    "Percent_With_At_Least_Six_Measurements"
]

print(
    cluster_summary[display_columns]
    .sort_values("All_Genes_Cluster_20")
    .to_string(index=False)
)


print("\n========== ALL-GENE CLUSTERING COMPLETE ==========")

print("\nImportant interpretation notes:")
print(
    "1. No master-matrix rows were excluded before clustering."
)
print(
    "2. Missing values were replaced by zero only in the temporary "
    "clustering representation."
)
print(
    "3. Zero means the average standardized response; it does not prove "
    "that an unmeasured gene had no biological response."
)
print(
    "4. Use Measured_Omics_Count and Information_Level when assessing "
    "the reliability of individual cluster assignments."
)
print(
    "5. The original normalized master matrix was not modified."
)