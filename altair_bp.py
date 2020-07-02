import pandas as pd
import altair as alt
import seaborn as sns
import streamlit as st
import numpy as np

st.title("Exploratory Dashboard")

@st.cache
def load_data_metrics_f1():
    df = pd.read_csv("data/hiv_protease/benchmark/metrics/f1.csv", index_col=0)
    df["fold"] = df.index
    return df

df = load_data_metrics_f1()

values = st.sidebar.slider('Select a range of values', 0.0, 1.0, (0.75, 1.0))
min, max = values

medians = df.apply(np.median)
indices = medians[(min < medians) & (medians < max)].index
medians = medians.to_frame(name="median_f1")
medians["Encoding"] = medians.index

medians = medians.drop("fold").loc[indices, :]

selectionc1 = alt.selection_interval(encodings=["x"])

c1 = alt.Chart(medians).mark_point().encode(
    x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
    y="median_f1",
    opacity=alt.condition(selectionc1, alt.value(1.0), alt.value(0.1)),
).add_selection(
    selectionc1
).properties(width=300)

st.sidebar.altair_chart(c1)


# print(df.head())

# print(indices)
df = df.loc[:, indices.tolist()+["fold"]]

print(df.shape)

@st.cache
def load_dfm(df):

    dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
                  var_name="Encoding", value_name="F1")

    def is_struc_based(e):
        if "asa" in e:
            return True
        elif "delaunay" in e:
            return True
        elif "disorder" in e:
            return True
        elif "hull" in e:
            return True
        elif "qsar" in e:
            return True
        elif "sse" in e:
            return True
        elif "ta" in e:
            return True
        else:
            return False

    dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["Encoding"]]

    return dfm

dfm = load_dfm(df)

c2 = alt.Chart(dfm).mark_boxplot().encode(
    y=alt.Y("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
    x=alt.X("F1"),
    # opacity=alt.condition(selectionc1, alt.value(1.0), alt.value(0.1)),
    color="type:N"
).transform_filter(
    selectionc1
).properties(width=600, height=800)

st.altair_chart(alt.hconcat(c1, c2))

# selection = alt.selection_interval(encodings=["x"])
# selection2 = alt.selection_interval(encodings=["x"])
# selection3 = alt.selection_multi(fields=['type'], bind='legend')
# selection4 = alt.selection_single()

# first = \
#     alt.Chart(dfm).mark_point().encode(
#         x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
#         y="average(F1)",
#         # opacity=alt.condition(selection, alt.value(1.0), alt.value(0.1)),
#         color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
#         tooltip="Encoding:N"
#     ).properties(width=600)
       # .add_selection(
       #  selection
    # )
        # .add_selection(
        # selection3
# )
#     .transform_filter(
#         selection3
#     ).properties(width=600)

# second = \
#     alt.Chart(dfm).mark_point().encode(
#         x=alt.X("Encoding", sort="-y"), #, axis=alt.Axis(labels=False, ticks=False, title=None)),
#         y="average(F1)",
#         color="type:N"
#     ).add_selection(
#         selection2
#     ).transform_filter(
#         selection
#     ).properties(width=600)
#
# third = \
#     alt.Chart(dfm).mark_boxplot().encode(
#         y=alt.Y("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
#         x=alt.X("F1"),
#         opacity=alt.condition(selection2, alt.value(1.0), alt.value(0.1)),
#         color="type:N"
#     ).transform_filter(
#         selection
#     ).add_selection(selection4).properties(width=600, height=800)

# st.altair_chart(alt.hconcat(alt.vconcat(first, second), third), use_container_width=True)
# st.altair_chart(alt.hconcat(alt.vconcat(first)))

# chart = alt.hconcat(alt.vconcat(
#     # alt.Chart(dfm).transform_aggregate(
#     #     mean_f1='average(F1)',
#     #     groupby=["Encoding"],
#     # ).mark_bar().add_selection(
#     #     selection
#     # ).encode(
#     #     x=alt.X("mean_f1:N", bin=alt.Bin(extent=[0.0, 1.0], step=0.01)),
#     #     y=alt.Y("count():Q")
#     # ).properties(width=1000),
#     # alt.Chart(dfm).mark_bar().add_selection(
#     #     selection
#     # ).encode(
#     #     x=alt.X("group:Q"),
#     #     y=alt.Y("count():Q")
#     # )
#     first, second,
#         alt.Chart(dfm).mark_bar().encode(
#             x=alt.X("F1", bin=alt.Bin(maxbins=20)),
#             y="count()",
#             opacity=alt.condition(selection4, alt.value(1.0), alt.value(0.1)),
#         ).transform_filter(selection4)
# ), third
# )

# df = pd.read_csv("data/hiv_protease/benchmark/similarity/seq_vs_str/diversity.csv", index_col=0)
#
# g2 = sns.clustermap(df.values)
# heatmap_data = g2.data2d
#
# import numpy as np
#
# x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
#
# source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "Diversity": heatmap_data.values.ravel()})
#
# e1, e2 = [], []
# for n, s in source.iterrows():
#     idx = heatmap_data.index[int(s["y"])]
#     col = heatmap_data.columns[int(s["x"])]
#     e1 += [df.index[idx]]
#     e2 += [df.columns[col]]
#
# source["Encoding1"] = e1
# source["Encoding2"] = e2
#
# # print(source.iloc[0, :])
# # print(heatmap_data.iloc[0,0])
#
# hm = alt.Chart(source).mark_rect().encode(
#     x=alt.X('x:O', title="Encoding 2", axis=alt.Axis(labels=False, ticks=False)),
#     y=alt.X('y:O', title="Encoding 1", axis=alt.Axis(labels=False, ticks=False)),
#     color='Diversity:Q',
#     tooltip=["Encoding1", "Encoding2", "Diversity"]
# ).properties(
#     title="Diversity",
#     width=600,
#     height=600
# )
#
# hm = hm.add_selection(alt.selection_single())
#
# chart = alt.vconcat(
#     alt.hconcat(alt.vconcat(first, second), third),
#     hm
# ).properties(
#     title="HIV protease"
# ).configure_title(
#     anchor="middle",
#     fontSize=26,
#     color="purple"
# )
#
# # import json
# # import pprint
# # js = json.loads(chart.to_json())
# # js["datasets"] = {}
# # pprint.pprint(js)
# chart.save("bp2.html")