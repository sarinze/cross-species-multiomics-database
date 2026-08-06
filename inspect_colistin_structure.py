import pandas as pd

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 1000)

df = pd.read_csv(
    "PA14_Colistin_Full_RNASeq_Data_Info.csv",
    header=None
)

print(df.iloc[0:10, :25])