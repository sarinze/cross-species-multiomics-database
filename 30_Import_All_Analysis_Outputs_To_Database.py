import sqlite3
import pandas as pd
import os

DATABASE = "pseudomonas_new.db"

FILES_TO_IMPORT = {
    # Step 23: Cluster validation
    "Cluster_Validation_Metrics.csv": "Cluster_Validation_Metrics",

    # Step 24/25: Final clustered matrix and heatmap outputs
    "Final_Clustered_Matrix.csv": "Final_Clustered_Matrix",
    "Cluster_Profile_Matrix_Mean.csv": "Cluster_Profile_Matrix_Mean",
    "Cluster_Heatmap_Data_Long.csv": "Cluster_Heatmap_Data_Long",

    # Step 26/27: PCA outputs
    "PCA_Outputs/PCA_Gene_Coordinates.csv": "PCA_Gene_Coordinates",
    "PCA_Outputs/PCA_Cluster_Centroids.csv": "PCA_Cluster_Centroids",
    "PCA_Outputs/PCA_Loadings.csv": "PCA_Loadings",
    "PCA_Outputs/PCA_Explained_Variance.csv": "PCA_Explained_Variance",

    # Step 28: Conserved response outputs
    "Conserved_Response_Analysis_Outputs/Conserved_Response_Genes.csv": "Conserved_Response_Genes",
    "Conserved_Response_Analysis_Outputs/Divergent_Response_Genes.csv": "Divergent_Response_Genes",
    "Conserved_Response_Analysis_Outputs/Top_Conserved_Response_Candidates.csv": "Top_Conserved_Response_Candidates",
    "Conserved_Response_Analysis_Outputs/Conserved_Response_By_Cluster.csv": "Conserved_Response_By_Cluster",
}

conn = sqlite3.connect(DATABASE)

print("\n========== IMPORTING ANALYSIS OUTPUTS ==========\n")

for filepath, table_name in FILES_TO_IMPORT.items():

    if not os.path.exists(filepath):
        print(f"SKIPPED - File not found: {filepath}")
        continue

    try:
        df = pd.read_csv(filepath)

        df.to_sql(
            table_name,
            conn,
            if_exists="replace",
            index=False
        )

        print(f"IMPORTED: {filepath}")
        print(f"Table name: {table_name}")
        print(f"Rows: {len(df)}")
        print(f"Columns: {len(df.columns)}\n")

    except Exception as e:
        print(f"FAILED: {filepath}")
        print(f"Reason: {e}\n")

conn.close()

print("========== IMPORT COMPLETE ==========")
print("All available analysis outputs have been imported into the database.")