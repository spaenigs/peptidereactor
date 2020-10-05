from more_itertools import chunked
from iFeature import AAC
from scipy.spatial import ConvexHull
from modlamp.core import read_fasta
from sklearn.manifold import TSNE

import pandas as pd
import altair as alt

import re
import secrets

from glob import glob

fastas = glob("data/*/seqs.fasta")
classes = glob("data/*/classes.txt")

df_res_1 = pd.DataFrame()
df_res_2 = pd.DataFrame()
for ds_path in glob("data/*/csv/"):
    try:
        df1 = pd.read_csv(ds_path + "/distance_frequency/dist_freq_dn_5_dc_50.csv", index_col=0)
        df1 = df1.loc[df1["y"] == 1, :].drop("y", axis=1)
        df2 = pd.read_csv(ds_path + "/cksaagp/cksaagp_gap_2.csv", index_col=0)
        df2 = df2.loc[df2["y"] == 1, :].drop("y", axis=1)
    except Exception:
        continue
        # df2 = pd.read_csv(ds_path + "/cksaagp/cksaagp_gap_1.csv", index_col=0)
        # df2 = df2.loc[df2["y"] == 1, :].drop("y", axis=1)
    df1["dataset"] = re.findall("data/(.*?)/", ds_path)[0]
    df_res_1 = pd.concat([df_res_1, df1])
    df2["dataset"] = re.findall("data/(.*?)/", ds_path)[0]
    df_res_2 = pd.concat([df_res_2, df2])

X_embedded_1 = TSNE(n_components=1, n_jobs=10).fit_transform(df_res_1.iloc[:, :-1].values)
X_embedded_2 = TSNE(n_components=1, n_jobs=10).fit_transform(df_res_2.iloc[:, :-1].values)

import os

path = "tsne_data.csv"

if os.path.exists(path):
    df_tsne = pd.read_csv(path, index_col=0)
else:
    df_tsne = pd.DataFrame({"x": X_embedded_1[:, 0], "y": X_embedded_2[:, 0]})
    df_tsne["dataset"] = df_res_1["dataset"].to_list()
    df_tsne.to_csv(path)

group = lambda enc: enc[:3]
df_tsne["group"] = [group(x) for x in df_tsne["dataset"]]

alt.Chart(df_tsne).mark_point(filled=True).encode(
    x="x:Q",
    y="y:Q",
    color="group:N",
    tooltip="dataset:N"
).save("chart.html")




#
# for fasta_path, class_path in zip(fastas, classes):
#     seqs, names = read_fasta(fasta_path)
#     with open(class_path) as f:
#         classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
#     seq_tuples = [[name, tup[0]]
#                   for name, tup in zip(names, zip(seqs, classes))
#                   if tup[1] == 1]
#     df_tmp = pd.DataFrame([res[1:] for res in AAC.AAC(seq_tuples, order=None)][1:])
#     df_tmp["dataset"] = re.findall("data/(.*?)/", fasta_path)[0]
#     df_res = pd.concat([df_res, df_tmp])
#
# # df_res.columns = [res[1:] for res in AAC.AAC(seq_tuples, order=None)][0] + [df_res.columns[-1]]
# X_embedded = TSNE(n_components=2, n_jobs=10).fit_transform(df_res.iloc[:, :-1].values)
#
# df_tsne = pd.DataFrame(X_embedded)
# df_tsne.columns = ["x", "y"]
# df_tsne["dataset"] = df_res["dataset"].to_list()