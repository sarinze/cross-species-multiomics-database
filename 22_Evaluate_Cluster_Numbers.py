import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

INPUT_MATRIX = "Clustering_Input_Matrix.csv"
NORMALIZED_COLUMNS_FILE = "Normalized_Columns.txt"

OUTPUT_FILE = "Cluster_Number_Evaluation.csv"

CLUSTER_NUMBERS = [4, 6, 8, 10, 12, 15, 20]

print("\n========== STEP 22: EVALUATE CLUSTER NUMBERS ==========")

df = pd.read_csv(INPUT_MATRIX, low_memory=False)

with open(NORMALIZED_COLUMNS_FILE, "r") as f:
    omics_cols = [line.strip() for line in f if line.strip()]

print(f"\nLoaded clustering matrix: {INPUT_MATRIX}")
print(f"Rows: {df.shape[0]}")
print(f"Omics columns used: {len(omics_cols)}")

X = df[omics_cols].copy()

print("\nCalculating hierarchical clustering once...")
distance_matrix = pdist(X, metric="euclidean")
linkage_matrix = linkage(distance_matrix, method="ward")

evaluation_rows = []

for n_clusters in CLUSTER_NUMBERS:
    labels = fcluster(
        linkage_matrix,
        t=n_clusters,
        criterion="maxclust"
    )

    temp = pd.DataFrame({
        "Cluster": labels
    })

    cluster_sizes = temp["Cluster"].value_counts().sort_index()

    largest_cluster = cluster_sizes.max()
    smallest_cluster = cluster_sizes.min()
    median_cluster = cluster_sizes.median()
    number_small_clusters = (cluster_sizes < 20).sum()

    evaluation_rows.append({
        "N_Clusters": n_clusters,
        "Largest_Cluster_Size": largest_cluster,
        "Smallest_Cluster_Size": smallest_cluster,
        "Median_Cluster_Size": median_cluster,
        "Small_Clusters_Under_20_Genes": number_small_clusters,
        "Largest_Cluster_Percent": round((largest_cluster / len(df)) * 100, 2)
    })

    print(f"\n===== {n_clusters} clusters =====")
    print(cluster_sizes)

evaluation = pd.DataFrame(evaluation_rows)
evaluation.to_csv(OUTPUT_FILE, index=False)

print(f"\nSaved cluster number evaluation: {OUTPUT_FILE}")

print("\nSummary:")
print(evaluation)

print("\n========== CLUSTER NUMBER EVALUATION COMPLETE ==========")