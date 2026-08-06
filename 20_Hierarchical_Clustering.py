import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt

INPUT_FILE = "Master_Multiomics_Matrix_Normalized.csv"
COLUMNS_FILE = "Normalized_Columns.txt"

OUTPUT_MATRIX = "Clustering_Input_Matrix.csv"
OUTPUT_GENE_LIST = "Clustering_Gene_List.csv"
OUTPUT_DENDROGRAM = "Hierarchical_Clustering_Dendrogram.png"

print("\n========== STEP 20: HIERARCHICAL CLUSTERING ==========")

df = pd.read_csv(INPUT_FILE, low_memory=False)

with open(COLUMNS_FILE, "r") as f:
    omics_cols = [line.strip() for line in f if line.strip()]

print(f"\nLoaded normalized matrix: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")
print(f"Omics columns used for clustering: {len(omics_cols)}")

# Keep only rows with at least 3 measured omics values
df["Measured_Omics_Count"] = df[omics_cols].notna().sum(axis=1)
cluster_df = df[df["Measured_Omics_Count"] >= 3].copy()

print(f"\nRows before filtering: {df.shape[0]}")
print(f"Rows after keeping genes with >=3 omics values: {cluster_df.shape[0]}")

# Create readable gene labels
cluster_df["Gene_Label"] = (
    cluster_df["Pseudomonas_Original_Locus_Tag"].fillna("").astype(str)
    + " | "
    + cluster_df["Pseudomonas_Gene_Name"].fillna("").astype(str)
    + " | "
    + cluster_df["Klebsiella_Locus_Tag"].fillna("").astype(str)
)

# Build clustering matrix
X = cluster_df[omics_cols].copy()

# For clustering only: fill missing values with 0
# Since data are z-scored, 0 means average response for that column
X_filled = X.fillna(0)

# Remove rows with no variation across selected features
row_variance = X_filled.var(axis=1)
cluster_df = cluster_df[row_variance > 0].copy()
X_filled = X_filled.loc[cluster_df.index]

print(f"Rows after removing zero-variance rows: {cluster_df.shape[0]}")

# Save clustering input
clustering_output = pd.concat(
    [
        cluster_df[
            [
                "Ortholog_Pair_ID",
                "Pseudomonas_Original_Locus_Tag",
                "Pseudomonas_Gene_Name",
                "Klebsiella_Locus_Tag",
                "Klebsiella_Gene_Name",
                "Measured_Omics_Count",
                "Gene_Label"
            ]
        ],
        X_filled
    ],
    axis=1
)

clustering_output.to_csv(OUTPUT_MATRIX, index=False)

cluster_df[
    [
        "Ortholog_Pair_ID",
        "Pseudomonas_Original_Locus_Tag",
        "Pseudomonas_Gene_Name",
        "Klebsiella_Locus_Tag",
        "Klebsiella_Gene_Name",
        "Measured_Omics_Count"
    ]
].to_csv(OUTPUT_GENE_LIST, index=False)

print(f"\nSaved clustering matrix: {OUTPUT_MATRIX}")
print(f"Saved clustering gene list: {OUTPUT_GENE_LIST}")

# Safety check
if X_filled.shape[0] < 2:
    raise ValueError("Not enough rows for hierarchical clustering after filtering.")

# Hierarchical clustering
distance_matrix = pdist(X_filled, metric="euclidean")
linkage_matrix = linkage(distance_matrix, method="ward")

# Plot dendrogram
plt.figure(figsize=(14, 8))

dendrogram(
    linkage_matrix,
    no_labels=True,
    color_threshold=None
)

plt.title("Hierarchical Clustering of Normalized Multi-Omics Gene Responses")
plt.xlabel("Genes / Ortholog Pairs")
plt.ylabel("Distance")

plt.tight_layout()
plt.savefig(OUTPUT_DENDROGRAM, dpi=300)
plt.close()

print(f"Saved dendrogram: {OUTPUT_DENDROGRAM}")

print("\n========== CLUSTERING COMPLETE ==========")
print("Important:")
print("The normalized master file was not changed.")
print("Missing values were filled with 0 only in the clustering input matrix.")
print("Only genes with at least 3 measured omics values were clustered.")