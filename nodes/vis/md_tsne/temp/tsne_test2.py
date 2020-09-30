from more_itertools import chunked
from iFeature import AAC
from scipy.spatial import ConvexHull
from modlamp.core import read_fasta
from glob import glob
from sklearn.manifold import TSNE

import pandas as pd
import altair as alt

import re

# paths = glob("data/*/csv/distance_frequency/dist_freq_dn_100_dc_100.csv")
#
# df_res = pd.DataFrame()
# for p in paths:
#     df_tmp = pd.read_csv(p, index_col=0)
#     df_tmp = df_tmp.loc[df_tmp["y"] == 1, :]
#     df_tmp.drop("y", axis=1, inplace=True)
#     df_tmp["dataset"] = re.findall("data/(.*?)/", p)[0]
#     df_res = pd.concat([df_res, df_tmp])
#
# print(df_res.shape)
#
# X_embedded = TSNE(n_components=2, n_jobs=10).fit_transform(df_res.iloc[:, :-1].values)
#
# df_tsne = pd.DataFrame(X_embedded)
# df_tsne.columns = ["x", "y"]
# df_tsne["dataset"] = df_res["dataset"].to_list()
#
# df_tsne.to_csv("datasets_tsne2.csv")
df_tsne = pd.read_csv("../../../../datasets_tsne2.csv", index_col=0)

group = lambda enc: \
             "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:3]

df_tsne["group"] = [group(x) for x in df_tsne["dataset"]]

alt.Chart(df_tsne).mark_point(filled=True, size=15).encode(
    x="x:Q",
    y="y:Q",
    color="group:N"
).save("chart_tsne2.html")