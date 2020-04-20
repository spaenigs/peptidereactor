from plotly.subplots import make_subplots

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

df_blast = pd.read_csv("blast_test.csv", sep=",", index_col=0)
# df_blast["lengths"] = [int(r.split("_")[1]) for r in df_blast.index]
# df_blast["lengths_2"] = range(0, len(df_blast.index))

# TODO find better way, since np.abs(np.log(0.000000000000000000000000000000000000000000000009)) == 108.33
df_blast["evalue_2"] = [0 if i == 0 else np.abs(np.log(i)) for i in df_blast.evalue]

print(df_blast.loc[df_blast.qlen == 8750, :].sort_values(by="evalue").iloc[0, -5:-3])

fig = go.Figure()

df_blast = df_blast.loc[df_blast.qlen > 5000, :]

df_blast.qlen = [str(i) for i in df_blast.qlen]




for i in set(df_blast.qlen):
    fig.add_trace(
        go.Box(y=df_blast.loc[df_blast.qlen == i, "evalue"], name=str(i))
    )



# # Create figure with secondary y-axis
# fig = make_subplots(specs=[[{"secondary_y": True}]])
#
# tl_type = "ols" # "ols"/"lowess"
# tmp = px.scatter(df_blast, x="length_2", y="evalues_2", trendline=tl_type, trendline_color_override="rgba(208,28,139,1)")
# scatter = tmp.data[0]
# scatter.marker.color = "rgba(208,28,139,0.1)"
# trendline = tmp.data[1]
#
# fig.add_trace(
#     scatter,
#     secondary_y=False,
# )

# fig.add_trace(trendline)

# indices, no_hits = [], []
# for i in df_blast.length:
#     row_idx = int(np.median(df_blast.filter(like=str(i), axis=0)["length_2"].values))
#     indices += [row_idx]
#     no_hits += [df_blast.iloc[row_idx, 1]]
#
# sorted_tuples = sorted(zip(indices, no_hits), key=lambda tup: tup[0])
# indices_sorted, no_hits_sorted = zip(*sorted_tuples)
#
# fig.add_trace(
#     go.Scatter(x=indices_sorted, y=no_hits_sorted, name="# of no hits", marker_color="rgba(77,172,38,1)"),
#     secondary_y=True,
# )
#
fig.update_layout(
    xaxis_title="Sequence length",
    # yaxis_title="log(e-value) * (-1): the smaller, the better",
    yaxis_type="log",
    font=dict(
        size=18,
    )
)

fig.show()