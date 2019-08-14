import json
import re
import pandas as pd
import numpy as np
import itertools

df = pd.read_csv(str(snakemake.input), index_col=0)

# TODO blosum62
nodes = [{"id": idx, "group": re.match("(.*?)(\d{1,2}[ABC]?)?encoder.*", idx).group(1)}
         for idx in df.index]

groups = set(list(map(lambda idx: re.match("(.*?)(type\d{1,2}[ABC]?)?encoder.*", idx).group(1), df.index)))

links = []
for g1, g2 in itertools.combinations(groups, 2):
    df_filtered = df.filter(regex=f"^{g1}(type\d{{1,2}}[ABC]?)?encoder.*", axis=0)  # filter by rows
    df_filtered = df_filtered.filter(regex=f"^{g2}(type\d{{1,2}}[ABC]?)?encoder.*")  # filter by columns
    min_index_idx, min_col_idx, val = \
        sorted([(n1, n2, df_filtered.loc[n1, n2])
                for n1, n2 in zip(df_filtered.index, df_filtered.idxmin(axis=1).values)],
               key=lambda triple: triple[2])[0]
    if val <= 0.05:
        links += [{"source": min_index_idx, "target": min_col_idx, "pvalue": val, "within": 0}]

for group in groups:
    df_f2 = df.filter(regex=f"^{group}(type\d{{1,2}}[ABC]?)?encoder.*", axis=0)
    df_f2 = df_f2.filter(regex=f"^{group}(type\d{{1,2}}[ABC]?)?encoder.*", axis=1)
    if np.sum(df_f2.shape) > 2:
        for n1, n2 in itertools.combinations(df_f2.index, 2):
            links += [{"source": n1, "target": n2, "pvalue": df_f2.loc[n1, n2], "within": 1}]

res = json.dumps({"nodes": nodes, "links": links})

with open(str(snakemake.output), mode="a") as f:
    f.write(res)
    f.flush()
