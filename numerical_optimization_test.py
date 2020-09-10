from sklearn.decomposition import PCA
from pathos.multiprocessing import ProcessingPool as Pool
from more_itertools import chunked
from functools import partial
from glob import glob

import pandas as pd
import altair as alt

import numpy as np
aaindex = pd.read_csv("peptidereactor/iFeature/data/AAindex.txt", sep="\t", index_col=0)
aaindex.columns = aaindex.columns[1:].tolist() + ["NaN"]
aaindex = aaindex.iloc[:, :-1].transpose()

df = aaindex.corr()

# referenz: Mittelpunkt, d.h. encoding bei ~0.3
# diversity: nach oben (Richtung 1.0) div steit, nach unten (Richtung 0.0) div sinkt
# phi corr: nach oben (Richtung 1.0) phi steigt, nach unten (Richtung 0.0) phi sinkt

df1 = pd.read_csv("data/ace_vaxinpad/benchmark/similarity/all_vs_all/diversity.csv", index_col=0)

pca = PCA(n_components=1)
X_new = pca.fit_transform(df1)

encoding_dist1 = pd.Series([x[0] for x in X_new], index=df1.columns)
encoding_dist1 = ((encoding_dist1 - encoding_dist1.min())/(encoding_dist1.max()-encoding_dist1.min()))
encoding_dist1.sort_values(ascending=False, inplace=True)


df_res1 = pd.DataFrame()
for i in encoding_dist1.index:
    df_tmp = pd.DataFrame({i: [df1.loc["flgc_aaindex_QIAN880102", i]]})
    df_res1 = pd.concat([df_res1, df_tmp.transpose()])

df_res1.columns = ["dist"]

df_res1["x"] = "flgc_aaindex_QIAN880102"
df_res1["y"] = df_res1.index

c1 = alt.Chart(df_res1).mark_point(filled=True).encode(
    x=alt.X("y:N", title="diversity", axis=alt.Axis(labels=False), sort=alt.Sort(alt.SortArray(encoding_dist1.index.values))),
    y=alt.Y("x:N", title=None),
    size="dist:Q",
    tooltip="y:N"
).properties(width=700)

###

df2 = pd.read_csv("data/ace_vaxinpad/benchmark/similarity/all_vs_all/phi.csv", index_col=0)

pca = PCA(n_components=1)
X_new = pca.fit_transform(df2)

encoding_dist2 = pd.Series([x[0] for x in X_new], index=df2.columns)
encoding_dist2 = ((encoding_dist2 - encoding_dist2.min())/(encoding_dist2.max()-encoding_dist2.min()))
encoding_dist2.sort_values(ascending=False, inplace=True)

df_res2 = pd.DataFrame()
for i in encoding_dist2.index:
    df_tmp = pd.DataFrame({i: [df2.loc["flgc_aaindex_QIAN880102", i]]})
    df_res2 = pd.concat([df_res2, df_tmp.transpose()])

df_res2.columns = ["dist"]

df_res2["x"] = "flgc_aaindex_QIAN880102"
df_res2["y"] = df_res2.index

c2 = alt.Chart(df_res2).mark_point(filled=True).encode(
    x=alt.X("y:N", title="phi", axis=alt.Axis(labels=False), sort=alt.Sort(alt.SortArray(encoding_dist2.index.values))),
    y=alt.Y("x:N", title=None),
    size="dist:Q",
    tooltip="y:N"
).properties(width=700)





(c1 & c2).save("chart.html")

