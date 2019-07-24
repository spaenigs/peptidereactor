import numpy as np
import pandas as pd
import os
import glob
from shutil import copyfile

print("hi from script")

df = pd.read_csv(str(snakemake.input), index_col=0)
df = df.applymap(lambda x: 1.0 if np.isnan(x) else x).transpose() * df.applymap(lambda x: 1.0 if np.isnan(x) else x)

tmp = df.copy()

# get high correlated indices and sum up occurrence
tmp["sum"] = tmp.applymap(lambda e: 1 if e >= 0.8 else 0).apply(lambda row: row.sum(), axis=1)

# sort by indices, which have many high correlated indices
tmp = tmp.sort_values(by="sum", ascending=False)
tmp["type"] = -1

res = []
for i in range(tmp.shape[0]):

    # assign a specific class to correlated indices
    df_filtered = tmp.loc[tmp["type"] == -1]

    # stop algorithm as soon as all correlated indices are found
    if df_filtered.shape[0] == 0:
        break

    # get the index names of high correlated indices for the current amino acid index
    hit_idx = tmp.loc[tmp.index[i], np.abs(tmp.loc[tmp.index[i]]) >= 0.8].index

    # only keep those names, which have not already assigned to a type, i.e., those indices,
    # which have not shown any correlation to a seen amino acid index (aaindex)
    new_idx = list(set(df_filtered.index).intersection(set(hit_idx)))

    # remove current aaindex, due to R(aaindex[i], aaindex[i]) == 1.0
    try:
        new_idx.remove(tmp.index[i])
    except ValueError:
        pass

    if len(new_idx) >= 1:
        # sort by correlation, descending
        t = tmp.loc[tmp.index[i], new_idx].sort_values(ascending=False)
        # keep highest correlated index
        res += [t.index[0]]

    # assign type to correlated indices
    tmp.loc[new_idx, "type"] = i

    # remove columns from original dataset to avoid adding of already seen indices
    try:
        tmp.drop(columns=new_idx + [tmp.index[i]], inplace=True)
    except KeyError:
        pass

with open(str(snakemake.output), mode="a") as f:
    for r in res:
        f.write(f"{r}\n")
        f.flush()
        src = glob.glob(f"data/in/csv/aaindex/*{r}.csv")[0]
        dst = f"data/out/aaindex/csv/final/normalized-{snakemake.wildcards.normalized}/{os.path.basename(src)}"
        copyfile(src, dst)
