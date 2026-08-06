import pandas as pd

print("\n========== STEP 23e: SELECT FINAL CLUSTERS ==========\n")

INPUT_ASSIGNMENTS = "Filtered_Cluster_Assignments.csv"
INPUT_MATRIX = "Clustering_Input_Matrix_Filtered.csv"

SELECTED_CLUSTER_COLUMN = "Cluster_20"

OUTPUT_FILE = "Final_Cluster_Assignments.csv"
OUTPUT_MATRIX_WITH_CLUSTERS = "Final_Clustered_Matrix.csv"
OUTPUT_SIZE_SUMMARY = "Final_Cluster_Size_Summary.csv"

assignments = pd.read_csv(INPUT_ASSIGNMENTS)
matrix = pd.read_csv(INPUT_MATRIX)

print(f"Loaded assignments: {INPUT_ASSIGNMENTS}")
print(f"Loaded matrix: {INPUT_MATRIX}")

if SELECTED_CLUSTER_COLUMN not in assignments.columns:
    raise ValueError(f"{SELECTED_CLUSTER_COLUMN} not found in assignment file.")

final_assignments = assignments.copy()
final_assignments["Final_Cluster"] = final_assignments[SELECTED_CLUSTER_COLUMN]

# Keep metadata columns + final cluster
metadata_cols = [
    col for col in final_assignments.columns
    if not col.startswith("Cluster_")
]

final_assignments = final_assignments[metadata_cols]

# Add final cluster to filtered matrix
matrix_with_clusters = matrix.copy()
matrix_with_clusters["Final_Cluster"] = final_assignments["Final_Cluster"]

# Cluster size summary
size_summary = (
    matrix_with_clusters["Final_Cluster"]
    .value_counts()
    .sort_index()
    .reset_index()
)

size_summary.columns = ["Final_Cluster", "Cluster_Size"]
size_summary["Cluster_Percent"] = (
    size_summary["Cluster_Size"] / len(matrix_with_clusters)
) * 100

# Save outputs
final_assignments.to_csv(OUTPUT_FILE, index=False)
matrix_with_clusters.to_csv(OUTPUT_MATRIX_WITH_CLUSTERS, index=False)
size_summary.to_csv(OUTPUT_SIZE_SUMMARY, index=False)

print(f"\nSelected cluster solution: {SELECTED_CLUSTER_COLUMN}")
print(f"Saved final assignments: {OUTPUT_FILE}")
print(f"Saved final clustered matrix: {OUTPUT_MATRIX_WITH_CLUSTERS}")
print(f"Saved cluster size summary: {OUTPUT_SIZE_SUMMARY}")

print("\nFinal cluster size summary:")
print(size_summary.to_string(index=False))

print("\n========== STEP 23e COMPLETE ==========\n")