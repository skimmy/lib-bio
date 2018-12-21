import sys
import numpy as np
import pandas as pd

n = 9
cols = [str(x) for x in range(n+1)]
header = ["x"] + cols + ["Ecc", "Alpha"]
df = pd.read_csv(sys.argv[1], header=None, names=header)
tmp = df[cols]*np.arange(n+1)
df['Ecc'] = tmp.sum(axis=1) / np.power(4,n)
df['Alpha'] = df['Ecc'] / n
df_sort = df.sort_values(["Alpha"])
df.to_csv("eccentricity_n"+str(n)+"_post.csv", index=False)
df_sort[["x", "Alpha"]].to_csv("eccentricity_n"+str(n)+"_reduced.csv", index=False)
