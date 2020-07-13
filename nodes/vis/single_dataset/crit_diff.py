from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

import pandas as pd
import altair as alt
import numpy as np


def dot_chart(df_f1, df_cd, cd, dataset):

    def cluster(values, axis):
        if axis == 1:
            values = values.T
        linkage = hierarchy.linkage(pdist(values), method="average", metric="euclidean")
        return hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)["leaves"]

    row_indices = cluster(df_cd.values, axis=0)
    col_indices = cluster(df_cd.values, axis=1)
    heatmap_data = df_cd.iloc[row_indices, col_indices]

    x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
    source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "cd": heatmap_data.values.ravel()})
    source["cd_cat"] = source["cd"].apply((lambda x: "critical different" if np.abs(x) > cd else "no difference"))

    e1, e2 = [], []
    for n, s in source.iterrows():
        e1 += [df_cd.index[int(s["y"])]]
        e2 += [df_cd.columns[int(s["x"])]]

    source["Encoding1"], source["Encoding2"] = e1, e2

    dfm_count = df_f1.apply(np.mean).groupby(
        by=lambda x: "psekraac" if "lambda-corr" in x or "g-gap" in x else x[:6]).count().to_frame("count")
    dfm_count["group"] = dfm_count.index

    for i in dfm_count.index:

        if i == "psekraac":
            indices_filtered = [ii for ii in df_cd.index if "g-gap" in ii or "lambda" in ii]
        else:
            indices_filtered = [ii for ii in df_cd.index if i == ii[:6]]

        df_cd_sub = df_cd.loc[indices_filtered, indices_filtered]
        df_le = df_cd_sub.le(-cd)
        df_ge = df_cd_sub.ge(cd)

        r = sum(df_le.values[np.triu_indices_from(df_le.values)]) + \
            sum(df_ge.values[np.triu_indices_from(df_ge.values)])

        dfm_count.loc[i, "cd_count"] = r
        dfm_count.loc[i, "cd_count_max"] = \
            len(indices_filtered) * (len(indices_filtered) - 1) / 2

    d, r = ["critical different", "no difference"], ["black", "gainsboro"]

    hm = alt.Chart(source).mark_rect().encode(
        x=alt.X('x:O', axis=alt.Axis(title="Encoding 1", labels=False, ticks=False)),
        y=alt.X('y:O', axis=alt.Axis(title="Encoding 2", labels=False, ticks=False)),
        color=alt.Color(
            "cd_cat:N",
            scale=alt.Scale(domain=d, range=r),
            legend=alt.Legend(title="CD (p<0.01)")
        ),
        tooltip=["Encoding1", "Encoding2", "cd"]
    ).properties(
        height=600, width=600
    )

    hm = alt.hconcat(
        hm,
        title=alt.TitleParams(text="All (left) vs. within group (right)", anchor="middle")
    )

    dfm_count = dfm_count.loc[dfm_count["count"] > 1, :]

    bars1 = alt.Chart(dfm_count).mark_bar(color=r[1]).encode(
        x='cd_count_max:Q',
        y="group:O",
        tooltip=["cd_count", "cd_count_max"]
    )

    bars2 = alt.Chart(dfm_count).mark_bar(color=r[0]).encode(
        x=alt.X('cd_count:Q', title="# critical different"),
        y=alt.Y("group:O", title="Encoding group"),
        tooltip=["cd_count", "cd_count_max"]
    )

    bars = (bars1 + bars2).properties(height=740)

    dots_data = df_cd.values[np.triu_indices_from(df_cd.values)]
    df_dots = pd.DataFrame({"x": dots_data, "y": [dataset] * len(dots_data)})

    dots = alt.Chart(df_dots).mark_circle(color=r[0]).encode(
        x=alt.X("binned_cd:Q", title="Critical difference (binned)"),
        y=alt.Y("y:N", axis=alt.Axis(title=None)),
        size=alt.Size("count(binned_cd):Q", legend=None),
        color=alt.condition(np.abs(alt.datum.binned_cd) >= cd,
                            alt.value(r[0]),
                            alt.value(r[1])),
        tooltip=["count(binned_cd):Q", "binned_cd:Q"]
    ).transform_bin(
        "binned_cd", "x", bin=alt.Bin(extent=[-300, 300], step=30)
    ).properties(
        width=600,
        height=100
    )

    return (hm & dots) | bars


