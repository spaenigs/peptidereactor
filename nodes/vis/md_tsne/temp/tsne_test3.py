from more_itertools import chunked
from iFeature import AAC
from scipy.spatial import ConvexHull
from modlamp.core import read_fasta
from glob import glob
from sklearn.manifold import TSNE

import pandas as pd
import altair as alt

import re

# def tsne(paths, coln):
#     df_res = pd.DataFrame()
#     for p in paths:
#         df_tmp = pd.read_csv(p, index_col=0)
#         df_tmp = df_tmp.loc[df_tmp["y"] == 1, :]
#         df_tmp.drop("y", axis=1, inplace=True)
#         df_tmp["dataset"] = re.findall("data/(.*?)/", p)[0]
#         df_res = pd.concat([df_res, df_tmp])
#     X_embedded = TSNE(n_components=1, n_jobs=10).fit_transform(df_res.iloc[:, :-1].values)
#     df_tsne = pd.DataFrame(X_embedded)
#     df_tsne.columns = [coln]
#     df_tsne[f"dataset_{coln}"] = df_res["dataset"].to_list()
#     return df_tsne
#
#
# paths_1 = glob("data/*/csv/distance_frequency/dist_freq_dn_100_dc_100.csv")
# df_tsne_1 = tsne(paths_1, "x")
#
# paths_2 = glob("data/*/csv/dde.csv")
# df_tsne_2 = tsne(paths_2, "y")
#
#
# df_tsne = pd.concat([df_tsne_1, df_tsne_2], axis=1)
# df_tsne.drop("dataset_x", axis=1, inplace=True)
# print(df_tsne.head())
# df_tsne.columns = ["x", "y", "dataset"]

# df_tsne.to_csv("datasets_tsne3.csv")
df_tsne = pd.read_csv("../../../../datasets_tsne3.csv", index_col=0)

alt.Chart(df_tsne).mark_point(filled=True, size=15).encode(
    x="x:Q",
    y="y:Q",
    color="dataset:N"
).save("chart_tsne3.html")