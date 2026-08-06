import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

INPUT_MATRIX = "Clustering_Input_Matrix.csv"
NORMALIZED_COLUMNS_FILE = "Normalized_Columns.txt"

OUTPUT_FILE = "Cluster_Validation_Metrics.csv"

CLUSTER_NUMBERS = [4, 6, 8, 10, 12, 15, 20]

print("\n========== STEP 23: CLUSTER VALIDATION ==========")

df = pd.read_csv(INPUT_MATRIX, low_memory=False)

with open(NORMALIZED_COLUMNS_FILE, "r") as f:
    omics_cols = [line.strip() for line in f if line.strip()]

print(f"\nLoaded clustering matrix: {INPUT_MATRIX}")
print(f"Rows: {df.shape[0]}")
print(f"Omics columns used: {len(omics_cols)}")

X = df[omics_cols].copy()

print("\nCalculating hierarchical clustering...")
distance_matrix = pdist(X, metric="euclidean")
linkage_matrix = linkage(distance_matrix, method="ward")

results = []

for n_clusters in CLUSTER_NUMBERS:
    print(f"\nValidating {n_clusters} clusters...")

    labels = fcluster(
        linkage_matrix,
        t=n_clusters,
        criterion="maxclust"
    )

    n_actual_clusters = len(set(labels))

    if n_actual_clusters < 2:
        print("Skipped: fewer than 2 clusters.")
        continue

    silhouette = silhouette_score(X, labels)
    calinski = calinski_harabasz_score(X, labels)
    davies = davies_bouldin_score(X, labels)

    cluster_sizes = pd.Series(labels).value_counts()

    largest_size = cluster_sizes.max()
    smallest_size = cluster_sizes.min()
    median_size = cluster_sizes.median()
    largest_percent = round((largest_size / len(df)) * 100, 2)

    if largest_percent > 75:
        balance_status = "Poor balance - one dominant cluster"
        recommended = "No"
        interpretation = "Rejected because most genes fall into one large cluster."
    elif smallest_size <= 1:
        balance_status = "Over-fragmented - single-gene cluster present"
        recommended = "No"
        interpretation = "Rejected because clustering is beginning to fragment into very small clusters."
    elif largest_percent <= 50:
        balance_status = "Acceptable balance"
        recommended = "Yes"
        interpretation = "Suitable for downstream biological interpretation."
    else:
        balance_status = "Moderate balance"
        recommended = "Possible"
        interpretation = "May be usable, but still contains a large dominant cluster."

    results.append({
        "N_Clusters_Requested": n_clusters,
        "N_Clusters_Actual": n_actual_clusters,
        "Silhouette_Score": silhouette,
        "Calinski_Harabasz_Index": calinski,
        "Davies_Bouldin_Index": davies,
        "Largest_Cluster_Size": largest_size,
        "Smallest_Cluster_Size": smallest_size,
        "Median_Cluster_Size": median_size,
        "Largest_Cluster_Percent": largest_percent,
        "Cluster_Balance_Status": balance_status,
        "Recommended_For_Downstream": recommended,
        "Interpretation": interpretation
    })

results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_FILE, index=False)

print(f"\nSaved validation metrics: {OUTPUT_FILE}")

print("\nValidation summary:")
pd.set_option("display.max_columns", None)
print(results_df.to_string(index=False))

print("\nHow to interpret:")
print("Silhouette Score: higher is better.")
print("Calinski-Harabasz Index: higher is better.")
print("Davies-Bouldin Index: lower is better.")
print("Largest Cluster Percent: lower is usually better for biological interpretability.")
print("Recommended_For_Downstream is based on balance and biological usefulness, not only internal metrics.")

print("\n========== CLUSTER VALIDATION COMPLETE ==========")