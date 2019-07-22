from scipy.spatial.distance import euclidean

import pandas as pd


def write_empty_file(path):
    with open(path, mode="w") as f:
        f.write("")
        f.flush()


df = pd.read_csv(str(input), index_col=0)
df["name"] = df.index

df_type_filtered = df.loc[df.name.apply(lambda x: snakemake.wildcards.type in x), :]

names = df_type_filtered.name
points = dict(zip(names, map(lambda x: dict(zip(names, map(lambda x: 0.0, range(len(names))))),
                             range(len(names)))))

res = {}
for i in names:
    for j in names:
        if i != j:
            if points[i][j] == 0.0:
                # https://en.wikipedia.org/wiki/Geometric_median
                points[i][j] = euclidean(df_type_filtered.loc[i, ["x1", "x2"]],
                                         df_type_filtered.loc[j, ["x1", "x2"]])
            else:
                continue
    res[i] = [sum(sorted(points[i].values()))]

res_df = pd.DataFrame(res).transpose()

if not res_df.empty:
    res_df.columns = ["distance_gm"]

    tmp2 = pd.concat([df_type_filtered, res_df], axis=1)
    tmp2["min_gm"] = False

    actual_total_min = 1e6
    actual_total_min_name = ""
    for row in tmp2.loc[tmp2.name.apply(lambda x: snakemake.wildcards.type in x), :].iterrows():
        if row[1]["distance_gm"] < actual_total_min:
           actual_total_min = row[1]["distance_gm"]
           actual_total_min_name = row[0]
    tmp2.loc[actual_total_min_name, "min_gm"] = True

    tmp2.to_csv(str(snakemake.output))

else:
    write_empty_file(str(snakemake.output))