from functools import reduce, partial

import pandas as pd
import altair as alt
import numpy as np

import re

from modlamp.core import read_fasta
from glob import glob

from iFeature import AAC

fastas = glob("data/*/seqs.fasta")
classes = glob("data/*/classes.txt")

df_res = pd.DataFrame()
for fasta_path, class_path in zip(fastas, classes):
    seqs, names = read_fasta(fasta_path)
    with open(class_path) as f:
        classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
    seq_tuples = [[name, tup[0]] for name, tup in zip(names, zip(seqs, classes)) if tup[1] == 1]
    df_tmp = pd.DataFrame([res[1:] for res in AAC.AAC(seq_tuples, order=None)][1:])
    df_tmp["dataset"] = re.findall("data/(.*?)/", fasta_path)[0]
    df_res = pd.concat([df_res, df_tmp])

df_res.columns = [res[1:] for res in AAC.AAC(seq_tuples, order=None)][0] + [df_res.columns[-1]]

from sklearn.manifold import TSNE
X_embedded = TSNE(n_components=2, n_jobs=10).fit_transform(df_res.iloc[:, :-1].values)

df_tsne = pd.DataFrame(X_embedded)
df_tsne["dataset"] = df_res["dataset"].to_list()
df_tsne.to_csv("datasets_tsne.csv")

# def compute_median_f1(path):
#     dataset = re.findall("data/(.*?)/", path)[0]
#     return pd\
#         .read_csv(path, index_col=0)\
#         .apply(np.median)\
#         .to_frame(dataset)
#
#
# paths = \
#     "data/hiv_protease/benchmark/metrics/f1.csv", \
#     "data/ace_vaxinpad/benchmark/metrics/f1.csv"
#
# df_res = pd.concat(
#     axis=1,
#     objs=reduce(
#         lambda dfs, p: dfs + [compute_median_f1(p)],
#         paths, [])).\
#     sort_index().\
#     fillna(0.0).\
#     T
#
# x, y = np.meshgrid(df_res.columns, df_res.index)
#
# source = pd.DataFrame({"Encoding": x.ravel(),
#                        "Dataset": y.ravel(),
#                        "F1": df_res.values.ravel()
#                        })
#
# l = np.hstack((np.array([df_res.index]).T, df_res.values))
# l = sorted(l, key=lambda x: [*x[1:]], reverse=True)
# sorted_datasets = np.array(l)[:, 0]
#
# alt.Chart(source).mark_rect().encode(
#     x='Encoding:O',
#     y=alt.Y(
#         'Dataset:N',
#         sort=alt.Sort(alt.SortArray(sorted_datasets))
#     ),
#     color='F1:Q',
#     tooltip=["Encoding", "Dataset", "F1"]
# ).properties(
#     width=1700
# ).save("chart.html")