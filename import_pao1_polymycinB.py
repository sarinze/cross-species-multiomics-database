import pandas as pd
import sqlite3

df = pd.read_csv("PAO1_PolymycinB_Rifampicin_Table_Content.csv")
conn = sqlite3.connect("pseudomonas_new.db")

df.to_sql("PAO1_PolymycinB_Rifampicin_Content", conn, if_exists="replace", index=False)

conn.close()
print("Content table imported.")