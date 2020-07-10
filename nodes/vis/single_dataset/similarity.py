from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

import os

import altair as alt
import pandas as pd
import numpy as np


def _get_clustering(df):

    def cluster(values, axis):
        if axis == 1:
            values = values.T
        linkage = hierarchy.linkage(pdist(values), method="average", metric="euclidean")
        return hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)["leaves"]

    row_indices = cluster(df.values, axis=0)
    col_indices = cluster(df.values, axis=1)
    heatmap_data = df.iloc[row_indices, col_indices]

    x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
    source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "Similarity": heatmap_data.values.ravel()})

    e1, e2 = [], []
    for n, s in source.iterrows():
        e1 += [df.index[int(s["y"])]]
        e2 += [df.columns[int(s["x"])]]

    source["Encoding1"], source["Encoding2"] = e1, e2
    source["Similarity_cat"] = source.Similarity.apply(
        lambda x: "<= 0.5" if x <= 0.5 else "<= 0.7" if x <= 0.7 else "<= 1.0")

    return source


def _hm(source, show_x_title=True, show_y_title=True):
    d, r = ["<= 0.5", "<= 0.7", "<= 1.0"], ["white", "#9ecae1", "#3182bd"]
    x_config = alt.Axis(labels=False, ticks=False)
    y_config = alt.Axis(labels=False, ticks=False)
    if not show_x_title:
        x_config = alt.Axis(labels=False, ticks=False, title=None)
    if not show_y_title:
        y_config = alt.Axis(labels=False, ticks=False, title=None)
    return alt.Chart(source).mark_rect().encode(
        x=alt.X('x:O', title="Encoding 2", axis=x_config),
        y=alt.Y('y:O', title="Encoding 1", axis=y_config),
        color=alt.Color(
            "Similarity_cat:N",
            scale=alt.Scale(domain=d, range=r),
            legend=alt.Legend(title="Similarity Range")
        ),
        tooltip=["Encoding1", "Encoding2", "Similarity"]
    ).properties(
        width=600,
        height=600
    )


def similarity_chart(df_div, df_phi, dataset):

    # TODO remove in production
    p_div = f"div_all_vs_all_{dataset}.csv"
    if os.path.exists(p_div):
        df_div_clustered = pd.read_csv(p_div, index_col=0)
    else:
        df_div_1 = pd.read_csv(f"data/{dataset}/benchmark/similarity/all_vs_all/diversity.csv", index_col=0)
        df_div_clustered = _get_clustering(df_div_1)
        df_div_clustered.to_csv(p_div)

    p_phi = f"phi_all_vs_all_{dataset}.csv"
    if os.path.exists(p_phi):
        df_phi_clustered = pd.read_csv(p_phi, index_col=0)
    else:
        df_phi_1 = pd.read_csv(f"data/{dataset}/benchmark/similarity/all_vs_all/phi.csv", index_col=0)
        df_phi_clustered = _get_clustering(df_phi_1)
        df_phi_clustered.to_csv(p_phi)

    hm1 = _hm(df_div_clustered, show_x_title=False)
    hm2 = _hm(df_phi_clustered, show_x_title=False, show_y_title=False)
    hm3 = _hm(_get_clustering(df_div))
    hm4 = _hm(_get_clustering(df_phi), show_y_title=False)

    t1 = alt.hconcat(hm1, hm2, title=alt.TitleParams(text="All vs. all", anchor="middle"))
    t2 = alt.hconcat(hm3, hm4, title=alt.TitleParams(text="Seq. vs. struc.", anchor="middle"))

    return (t1 & t2).configure_view(strokeWidth=0)