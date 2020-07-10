from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

import pandas as pd
import altair as alt
import numpy as np


def dot_chart(df_f1, df_cd, cd):

    medians = df_f1.apply(np.median).sort_values(ascending=False)
    indices = medians.index.tolist()

    indices_filtered = [i for i in indices if "g-gap" in i or "lambda" in i]
    print(len(indices_filtered))
    # df_cd = df_cd.loc[indices_filtered, indices_filtered]

    print(df_cd.shape)
    # df_cd = df_cd.loc[indices[:100], indices[:100]]
    # print(df_cd.shape)

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

    input_dropdown = alt.binding_select(options=dfm_count.index.tolist())
    selection = alt.selection_single(name="e", fields=['Encoding1'], bind=input_dropdown)

    hm = alt.Chart(source).mark_rect().encode(
        x=alt.X('x:O', axis=alt.Axis(labels=False, ticks=False)),
        y=alt.X('y:O', axis=alt.Axis(labels=False, ticks=False)),
        color=alt.Color(
            "cd_cat:N",
            scale=alt.Scale(domain=["critical different", "no difference"], range=["black", "whitesmoke"]),
            legend=alt.Legend(title="CD (p<0.01)")
        ),
        tooltip=["Encoding1", "Encoding2", "cd"]
    ).properties(
        height=600, width=600
    )

    # hm = hm.add_selection(
    #     selection
    # )
    #
    # hm = hm.transform_filter(
    #     "(substring(datum.Encoding1, 0, 6) === e) && (substring(datum.Encoding2, 0, 6) === e)"
    # )

    return hm


