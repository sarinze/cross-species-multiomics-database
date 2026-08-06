import pandas as pd

pa14 = pd.read_csv("PA14_Genome_Annotation.csv")
mgh = pd.read_csv("Kleb_MGH_78578_Genome_Annotation.csv")

print("PA14")
print(pa14.columns.tolist())

print("\nMGH78578")
print(mgh.columns.tolist())