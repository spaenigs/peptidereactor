from plotly.subplots import make_subplots

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

df_blast = pd.read_csv("blast_test.csv", index_col=0)
df_blast["lengths"] = [int(r.split("_")[1]) for r in df_blast.index]
df_blast["lengths_2"] = range(0, len(df_blast.index))
df_blast["evalues_2"] = [0 if i == 0 else np.log(i) * (-1) for i in df_blast.evalues]

# Create figure with secondary y-axis
fig = make_subplots(specs=[[{"secondary_y": True}]])

tl_type = "lowess" # "ols"
tmp = px.scatter(df_blast, x="lengths_2", y="evalues_2", trendline=tl_type, trendline_color_override="rgba(208,28,139,1)")
scatter = tmp.data[0]
scatter.marker.color = "rgba(208,28,139,0.1)"
trendline = tmp.data[1]

fig.add_trace(
    scatter,
    secondary_y=False,
)

fig.add_trace(trendline)

indices, no_hits = [], []
for i in set([r.split("_")[1] for r in df_blast.index]):
    row_idx = int(np.median(df_blast.filter(like=f"len_{i}", axis=0)["lengths_2"].values))
    indices += [row_idx]
    no_hits += [df_blast.iloc[row_idx, 1]]

sorted_tuples = sorted(zip(indices, no_hits), key=lambda tup: tup[0])
indices_sorted, no_hits_sorted = zip(*sorted_tuples)

fig.add_trace(
    go.Scatter(x=indices_sorted, y=no_hits_sorted, name="# of no hits", marker_color="rgba(77,172,38,1)"),
    secondary_y=True,
)

fig.update_layout(
    xaxis_title="Sequence length",
    yaxis_title="log(e-value) * (-1): the smaller, the better",
    font=dict(
        size=18,
    )
)

fig.show()