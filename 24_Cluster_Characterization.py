import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist

print("\n========== STEP 24: CLUSTER CHARACTERIZATION ==========\n")

INPUT_FILE = "Final_Clustered_Matrix.csv"

OUTPUT_CLUSTER_SUMMARY = "Cluster_Characterization_Summary.csv"
OUTPUT_CLUSTER_PROFILES = "Cluster_Profile_Matrix.csv"
OUTPUT_GENE_LEVEL = "Cluster_Gene_Level_Characterization.csv"
OUTPUT_REPRESENTATIVE_GENES = "Cluster_Representative_Genes.csv"
OUTPUT_TOP_FEATURES = "Cluster_Top_Features.csv"
OUTPUT_EXCEL = "Cluster_Characterization_Workbook.xlsx"

df = pd.read_csv(INPUT_FILE)

print(f"Loaded: {INPUT_FILE}")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

if "Final_Cluster" not in df.columns:
    raise ValueError("Final_Cluster column not found. Run 23e first.")

metadata_keywords = [
    "gene", "locus", "tag", "name", "description",
    "product", "protein", "accession", "cluster"
]

metadata_cols = [
    col for col in df.columns
    if any(keyword in col.lower() for keyword in metadata_keywords)
]

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

omics_cols = [
    col for col in numeric_cols
    if col not in metadata_cols
    and col not in ["Final_Cluster", "Measured_Omics_Count"]
]

print(f"\nMetadata columns detected: {len(metadata_cols)}")
print(f"Omics columns used for characterization: {len(omics_cols)}")

# -----------------------------
# Cluster size annotation
# -----------------------------
def classify_cluster_size(size):
    if size >= 100:
        return "Major biological cluster"
    elif size >= 20:
        return "Minor biological cluster"
    elif size >= 5:
        return "Small candidate cluster"
    else:
        return "Outlier cluster"

cluster_sizes = (
    df["Final_Cluster"]
    .value_counts()
    .sort_index()
    .reset_index()
)

cluster_sizes.columns = ["Final_Cluster", "Cluster_Size"]
cluster_sizes["Cluster_Percent"] = (
    cluster_sizes["Cluster_Size"] / len(df)
) * 100

cluster_sizes["Cluster_Category"] = cluster_sizes["Cluster_Size"].apply(classify_cluster_size)

# -----------------------------
# Cluster profile matrix
# -----------------------------
profile_mean = df.groupby("Final_Cluster")[omics_cols].mean()
profile_median = df.groupby("Final_Cluster")[omics_cols].median()
profile_std = df.groupby("Final_Cluster")[omics_cols].std()

profile_mean.to_csv(OUTPUT_CLUSTER_PROFILES)

# -----------------------------
# Summary table
# -----------------------------
summary = cluster_sizes.copy()

for cluster_id in summary["Final_Cluster"]:
    cluster_data = df[df["Final_Cluster"] == cluster_id]

    values = cluster_data[omics_cols].values.flatten()
    values = values[~pd.isna(values)]

    summary.loc[summary["Final_Cluster"] == cluster_id, "Mean_All_Omics"] = np.mean(values)
    summary.loc[summary["Final_Cluster"] == cluster_id, "Median_All_Omics"] = np.median(values)
    summary.loc[summary["Final_Cluster"] == cluster_id, "Std_All_Omics"] = np.std(values)
    summary.loc[summary["Final_Cluster"] == cluster_id, "Min_All_Omics"] = np.min(values)
    summary.loc[summary["Final_Cluster"] == cluster_id, "Max_All_Omics"] = np.max(values)

summary.to_csv(OUTPUT_CLUSTER_SUMMARY, index=False)

# -----------------------------
# Gene-level characterization
# -----------------------------
gene_level = df.copy()

gene_level["Mean_Omics_Response"] = gene_level[omics_cols].mean(axis=1)
gene_level["Median_Omics_Response"] = gene_level[omics_cols].median(axis=1)
gene_level["Max_Omics_Response"] = gene_level[omics_cols].max(axis=1)
gene_level["Min_Omics_Response"] = gene_level[omics_cols].min(axis=1)
gene_level["Omics_Response_Range"] = (
    gene_level["Max_Omics_Response"] - gene_level["Min_Omics_Response"]
)

gene_level.to_csv(OUTPUT_GENE_LEVEL, index=False)

# -----------------------------
# Representative genes
# Genes closest to cluster centroid
# -----------------------------
X = df[omics_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
X_scaled = StandardScaler().fit_transform(X)

scaled_df = pd.DataFrame(X_scaled, columns=omics_cols)
scaled_df["Final_Cluster"] = df["Final_Cluster"].values

rep_rows = []

for cluster_id in sorted(df["Final_Cluster"].unique()):
    idx = df[df["Final_Cluster"] == cluster_id].index
    cluster_scaled = scaled_df.loc[idx, omics_cols]

    centroid = cluster_scaled.mean(axis=0).values.reshape(1, -1)
    distances = cdist(cluster_scaled.values, centroid).flatten()

    temp = df.loc[idx].copy()
    temp["Distance_To_Cluster_Centroid"] = distances

    temp = temp.sort_values("Distance_To_Cluster_Centroid").head(10)

    keep_cols = metadata_cols + [
        "Final_Cluster",
        "Distance_To_Cluster_Centroid",
        "Measured_Omics_Count"
    ]

    keep_cols = [col for col in keep_cols if col in temp.columns]

    rep_rows.append(temp[keep_cols])

representative_genes = pd.concat(rep_rows, ignore_index=True)
representative_genes.to_csv(OUTPUT_REPRESENTATIVE_GENES, index=False)

# -----------------------------
# Top positive and negative features per cluster
# -----------------------------
top_feature_rows = []

for cluster_id in profile_mean.index:
    means = profile_mean.loc[cluster_id].sort_values(ascending=False)

    top_positive = means.head(5)
    top_negative = means.tail(5)

    for feature, value in top_positive.items():
        top_feature_rows.append({
            "Final_Cluster": cluster_id,
            "Direction": "Highest mean feature",
            "Feature": feature,
            "Mean_Value": value
        })

    for feature, value in top_negative.items():
        top_feature_rows.append({
            "Final_Cluster": cluster_id,
            "Direction": "Lowest mean feature",
            "Feature": feature,
            "Mean_Value": value
        })

top_features = pd.DataFrame(top_feature_rows)
top_features.to_csv(OUTPUT_TOP_FEATURES, index=False)

# -----------------------------
# Excel workbook
# -----------------------------
with pd.ExcelWriter(OUTPUT_EXCEL, engine="openpyxl") as writer:
    summary.to_excel(writer, sheet_name="Cluster Summary", index=False)
    profile_mean.to_excel(writer, sheet_name="Mean Profiles")
    profile_median.to_excel(writer, sheet_name="Median Profiles")
    profile_std.to_excel(writer, sheet_name="Std Profiles")
    representative_genes.to_excel(writer, sheet_name="Representative Genes", index=False)
    top_features.to_excel(writer, sheet_name="Top Features", index=False)
    gene_level.to_excel(writer, sheet_name="Gene Level", index=False)

print(f"\nSaved: {OUTPUT_CLUSTER_SUMMARY}")
print(f"Saved: {OUTPUT_CLUSTER_PROFILES}")
print(f"Saved: {OUTPUT_GENE_LEVEL}")
print(f"Saved: {OUTPUT_REPRESENTATIVE_GENES}")
print(f"Saved: {OUTPUT_TOP_FEATURES}")
print(f"Saved: {OUTPUT_EXCEL}")

print("\nCluster summary:")
print(summary.to_string(index=False))

print("\n========== STEP 24 COMPLETE ==========\n")