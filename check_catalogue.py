import pandas as pd

file_name = "Final_Cluster_Gene_Catalogue.csv"

df = pd.read_csv(file_name, low_memory=False)

print("\nColumns in Final_Cluster_Gene_Catalogue.csv:\n")

for column in df.columns:
    print(column)

print("\nFirst five rows:\n")
print(df.head())