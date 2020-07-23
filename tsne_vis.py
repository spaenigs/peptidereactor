from sklearn.metrics import davies_bouldin_score

import pandas as pd
import altair as alt
import numpy as np

df_tsne = pd.read_csv("datasets_tsne.csv", index_col=0)
df_tsne["group"] = df_tsne["dataset"].apply(lambda e: e[:3])

# df_tsne.groupby().apply(lambda x: davies_bouldin_score(X, labels)

from scipy.spatial import ConvexHull, convex_hull_plot_2d
points = df_tsne.loc[df_tsne["dataset"].isin(["amp_csamp"]), :].iloc[:, :-2].values
hull = ConvexHull(points)
print(hull.area)
# print(df_tsne.head())

base = alt.Chart(df_tsne.loc[df_tsne["dataset"].isin(["amp_csamp"]), :])\

scatter = base.mark_circle().encode(
    x="0:Q",
    y="1:Q",
    color="group:N",
    # facet=alt.Facet('dataset:N', columns=6)
).resolve_scale(
    x="independent",
    y="independent"
).properties(
    height=100,
    width=100
)

hullc = alt.Chart().mark_point(color="red").encode(
    x="0:Q",
    y="1:Q"
)

# import matplotlib.pyplot as plt
# # plt.plot(points[:,0], points[:,1], 'o')
# for simplex in hull.simplices:
#     print((points[simplex, 0], points[simplex, 1]))
#
# print(hull.simplices)
#
# plt.show()

d = df_tsne.loc[df_tsne["dataset"].isin(["amp_csamp"]), :].iloc[hull.vertices, :-2]
source1 = pd.DataFrame({"x": d["0"], "y": d["1"]})
source1["order"] = range(1, source1.shape[0]+1)
source1.index = range(0, source1.shape[0])

print(source1)

hullc = alt.Chart(source1).mark_line().encode(
    x="x:Q",
    y="y:Q",
    order="order",
    tooltip=["order", "x", "y"]
)

# hullc.save("chart.html")

source = pd.DataFrame({
    "x": source1["x"].to_list() + [source1["x"].to_list()[0]],
    "y": source1["y"].to_list() + [source1["y"].to_list()[0]],
    "order": source1["order"].to_list() + [10]}
)

print(source)

c2 = alt.Chart(source).mark_line().encode(
    x="x:Q",
    y="y:Q",
    order="order"
)
(hullc | c2).save("chart.html")