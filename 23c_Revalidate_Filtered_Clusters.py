import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from scipy.cluster.hierarchy import linkage, fcluster

print("\n========== STEP 23c: REVALIDATE FILTERED CLUSTERS ==========\n")

INPUT_FILE = "Clustering_Input_Matrix_Filtered.csv"
OUTPUT_METRICS_FILE = "Filtered_Cluster_Validation_Metrics.csv"
OUTPUT_ASSIGNMENTS_FILE = "Filtered_Cluster_Assignments.csv"

CLUSTER_NUMBERS = [4, 6, 8, 10, 12, 15, 20]

# -----------------------------
# Load filtered matrix
# -----------------------------
df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# -----------------------------
# Identify metadata and omics columns
# -----------------------------
metadata_keywords = [
    "gene",
    "locus",
    "tag",
    "name",
    "description",
    "product",
    "protein",
    "accession",
    "cluster",
]

metadata_cols = [
    col for col in df.columns
    if any(keyword in col.lower() for keyword in metadata_keywords)
]

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

omics_cols = [
    col for col in numeric_cols
    if col not in metadata_cols
    and col != "Measured_Omics_Count"
]

print(f"\nMetadata columns detected: {len(metadata_cols)}")
print(f"Omics columns used: {len(omics_cols)}")

print("\nOmics columns:")
for col in omics_cols:
    print(f" - {col}")

# -----------------------------
# Prepare data
# -----------------------------
X = df[omics_cols].copy()
X = X.replace([np.inf, -np.inf], np.nan)
X = X.fillna(0)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# -----------------------------
# Hierarchical clustering
# -----------------------------
print("\nCalculating hierarchical clustering...")

Z = linkage(X_scaled, method="ward")

# -----------------------------
# Validate cluster numbers
# -----------------------------
metrics = []
assignment_df = df[metadata_cols].copy()

for n_clusters in CLUSTER_NUMBERS:
    print(f"\nValidating {n_clusters} clusters...")

    labels = fcluster(Z, t=n_clusters, criterion="maxclust")

    n_actual = len(np.unique(labels))

    silhouette = silhouette_score(X_scaled, labels)
    calinski = calinski_harabasz_score(X_scaled, labels)
    davies = davies_bouldin_score(X_scaled, labels)

    cluster_sizes = pd.Series(labels).value_counts().sort_index()

    largest = cluster_sizes.max()
    smallest = cluster_sizes.min()
    median = cluster_sizes.median()
    largest_percent = (largest / len(labels)) * 100

    if largest_percent > 80:
        balance_status = "Poor balance: one cluster dominates"
        recommended = "No"
    elif largest_percent > 60:
        balance_status = "Moderate imbalance"
        recommended = "Maybe"
    else:
        balance_status = "Reasonably balanced"
        recommended = "Yes"

    metrics.append({
        "N_Clusters_Requested": n_clusters,
        "N_Clusters_Actual": n_actual,
        "Silhouette_Score": silhouette,
        "Calinski_Harabasz_Index": calinski,
        "Davies_Bouldin_Index": davies,
        "Largest_Cluster_Size": largest,
        "Smallest_Cluster_Size": smallest,
        "Median_Cluster_Size": median,
        "Largest_Cluster_Percent": largest_percent,
        "Cluster_Balance_Status": balance_status,
        "Recommended_For_Downstream": recommended,
    })

    assignment_df[f"Cluster_{n_clusters}"] = labels

metrics_df = pd.DataFrame(metrics)

metrics_df.to_csv(OUTPUT_METRICS_FILE, index=False)
assignment_df.to_csv(OUTPUT_ASSIGNMENTS_FILE, index=False)

print(f"\nSaved validation metrics: {OUTPUT_METRICS_FILE}")
print(f"Saved cluster assignments: {OUTPUT_ASSIGNMENTS_FILE}")

print("\nValidation summary:")
print(metrics_df.to_string(index=False))

print("\n========== STEP 23c COMPLETE ==========\n")