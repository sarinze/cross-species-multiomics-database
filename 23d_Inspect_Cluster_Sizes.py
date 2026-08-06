import pandas as pd
import matplotlib.pyplot as plt

print("\n========== STEP 23d: INSPECT CLUSTER SIZES ==========\n")

INPUT_FILE = "Filtered_Cluster_Assignments.csv"
OUTPUT_SUMMARY = "Filtered_Cluster_Size_Inspection_Summary.csv"
OUTPUT_DETAIL = "Filtered_Cluster_Size_Details.csv"

df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

cluster_cols = [col for col in df.columns if col.startswith("Cluster_")]

print(f"\nCluster solutions found: {cluster_cols}")

summary_rows = []
detail_rows = []

for col in cluster_cols:
    n_requested = col.replace("Cluster_", "")

    cluster_sizes = df[col].value_counts().sort_index()

    total_genes = len(df)
    largest = cluster_sizes.max()
    smallest = cluster_sizes.min()
    median = cluster_sizes.median()
    mean = cluster_sizes.mean()

    singleton_clusters = (cluster_sizes == 1).sum()
    clusters_lt_5 = (cluster_sizes < 5).sum()
    clusters_lt_10 = (cluster_sizes < 10).sum()

    genes_in_singletons = cluster_sizes[cluster_sizes == 1].sum()
    genes_in_lt_5 = cluster_sizes[cluster_sizes < 5].sum()
    genes_in_lt_10 = cluster_sizes[cluster_sizes < 10].sum()

    largest_percent = (largest / total_genes) * 100
    singleton_gene_percent = (genes_in_singletons / total_genes) * 100
    lt5_gene_percent = (genes_in_lt_5 / total_genes) * 100
    lt10_gene_percent = (genes_in_lt_10 / total_genes) * 100

    summary_rows.append({
        "Cluster_Solution": col,
        "N_Clusters": len(cluster_sizes),
        "Total_Genes": total_genes,
        "Largest_Cluster_Size": largest,
        "Largest_Cluster_Percent": largest_percent,
        "Smallest_Cluster_Size": smallest,
        "Median_Cluster_Size": median,
        "Mean_Cluster_Size": mean,
        "Singleton_Clusters": singleton_clusters,
        "Clusters_With_Fewer_Than_5_Genes": clusters_lt_5,
        "Clusters_With_Fewer_Than_10_Genes": clusters_lt_10,
        "Genes_In_Singleton_Clusters": genes_in_singletons,
        "Genes_In_Clusters_Fewer_Than_5": genes_in_lt_5,
        "Genes_In_Clusters_Fewer_Than_10": genes_in_lt_10,
        "Percent_Genes_In_Singletons": singleton_gene_percent,
        "Percent_Genes_In_Clusters_Fewer_Than_5": lt5_gene_percent,
        "Percent_Genes_In_Clusters_Fewer_Than_10": lt10_gene_percent,
    })

    for cluster_id, size in cluster_sizes.items():
        detail_rows.append({
            "Cluster_Solution": col,
            "Cluster_ID": cluster_id,
            "Cluster_Size": size,
            "Cluster_Percent": (size / total_genes) * 100
        })

    plt.figure(figsize=(8, 5))
    plt.bar(cluster_sizes.index.astype(str), cluster_sizes.values)
    plt.xlabel("Cluster ID")
    plt.ylabel("Number of genes")
    plt.title(f"Cluster Size Distribution: {col}")
    plt.xticks(rotation=45)
    plt.tight_layout()

    plot_file = f"{col}_Size_Distribution.png"
    plt.savefig(plot_file, dpi=300)
    plt.close()

    print(f"Saved plot: {plot_file}")

summary_df = pd.DataFrame(summary_rows)
detail_df = pd.DataFrame(detail_rows)

summary_df.to_csv(OUTPUT_SUMMARY, index=False)
detail_df.to_csv(OUTPUT_DETAIL, index=False)

print(f"\nSaved summary: {OUTPUT_SUMMARY}")
print(f"Saved details: {OUTPUT_DETAIL}")

print("\nCluster size inspection summary:")
print(summary_df.to_string(index=False))

print("\n========== STEP 23d COMPLETE ==========\n")