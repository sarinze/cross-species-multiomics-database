import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

INPUT_MATRIX = "Clustering_Input_Matrix.csv"
NORMALIZED_COLUMNS_FILE = "Normalized_Columns.txt"

OUTPUT_FILE = "Hierarchical_Cluster_Assignments.csv"
SUMMARY_FILE = "Hierarchical_Cluster_Summary.csv"

N_CLUSTERS = 15

print("\n========== STEP 21: ASSIGN HIERARCHICAL CLUSTERS ==========")

df = pd.read_csv(INPUT_MATRIX, low_memory=False)

with open(NORMALIZED_COLUMNS_FILE, "r") as f:
    omics_cols = [line.strip() for line in f if line.strip()]

print(f"\nLoaded clustering matrix: {INPUT_MATRIX}")
print(f"Rows: {df.shape[0]}")
print(f"Omics columns used: {len(omics_cols)}")
print(f"Number of clusters requested: {N_CLUSTERS}")

X = df[omics_cols].copy()

distance_matrix = pdist(X, metric="euclidean")
linkage_matrix = linkage(distance_matrix, method="ward")

df["Hierarchical_Cluster"] = fcluster(
    linkage_matrix,
    t=N_CLUSTERS,
    criterion="maxclust"
)

# Save full assignment table
df.to_csv(OUTPUT_FILE, index=False)

# Create cluster summary
summary = (
    df.groupby("Hierarchical_Cluster")
    .agg(
        Gene_Count=("Ortholog_Pair_ID", "count"),
        Mean_Measured_Omics_Count=("Measured_Omics_Count", "mean")
    )
    .reset_index()
)

# Add average response per omics column for each cluster
cluster_means = df.groupby("Hierarchical_Cluster")[omics_cols].mean().reset_index()

summary = summary.merge(cluster_means, on="Hierarchical_Cluster", how="left")

summary.to_csv(SUMMARY_FILE, index=False)

print(f"\nSaved cluster assignments: {OUTPUT_FILE}")
print(f"Saved cluster summary: {SUMMARY_FILE}")

print("\nCluster sizes:")
print(summary[["Hierarchical_Cluster", "Gene_Count", "Mean_Measured_Omics_Count"]])

print("\n========== CLUSTER ASSIGNMENT COMPLETE ==========")