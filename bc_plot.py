from scipy import stats

import pandas as pd
import altair as alt
import numpy as np

import re

# dataset = "bce_bcm"
dataset = "hiv_protease"
# dataset = "ace_vaxinpad"

# df = pd.read_csv(f"bc_data_{dataset}.csv")
df1 = pd.read_csv(f"/home/spaenigs/Downloads/bc_data_hiv_protease.csv")
df1["dataset"] = "hiv_protease"
df2 = pd.read_csv(f"/home/spaenigs/Downloads/bc_data_ace_vaxinpad.csv")
df2["dataset"] = "ace_vaxinpad"
df3 = pd.read_csv(f"bc_data_bce_bcm.csv")
df3["dataset"] = "bachem (wl:7, co:1.5%)"

df = pd.concat([df1, df2, df3])

# for encoding_pair in df["group"].unique():
#     e1, e2 = re.findall("'(.*?)', '(.*?)'", encoding_pair)[0]
#     e1_resampled_f1s = df.loc[df["encoding"].str.contains("resampled"), :].loc[df["encoding"].str.contains(e1), "f1"].values
#     e1_non_resampled_f1s = df.loc[~df["encoding"].str.contains("resampled"), :].loc[df["encoding"].str.contains(e1), "f1"].values
#     _, e1_resampled_vs_non_pv = stats.ttest_rel(e1_non_resampled_f1s, e1_resampled_f1s)
#     print(e1_resampled_vs_non_pv)

bp1 = alt.Chart(df2).mark_boxplot().encode(
    x=alt.X("encoding:N", title=None, sort=df2.groupby("encoding").median().sort_values("f1", ascending=True).index.to_list()),
    y=alt.Y('f1:Q', scale=alt.Scale(domain=[0.0, 1.0])),
    color="resampled:N",
    row=alt.Row("dataset:N", title=None),
    column=alt.Column("group:N", header=alt.Header(labels=False))
).resolve_scale(
    x="independent",
).properties(
    width=120
)

bp2 = alt.Chart(df1).mark_boxplot().encode(
    x=alt.X("encoding:N", title=None, sort=df1.groupby("encoding").median().sort_values("f1", ascending=True).index.to_list()),
    y=alt.Y('f1:Q', scale=alt.Scale(domain=[0.0, 1.0])),
    color="resampled:N",
    row=alt.Row("dataset:N", title=None),
    column=alt.Column("group:N", header=alt.Header(labels=False))
).resolve_scale(
    x="independent",
)

bp3 = alt.Chart(df3).mark_boxplot().encode(
    x=alt.X("encoding:N", title=None, sort=df3.groupby("encoding").median().sort_values("f1", ascending=True).index.to_list()),
    y=alt.Y('f1:Q', scale=alt.Scale(domain=[0.0, 1.0])),
    color="resampled:N",
    row=alt.Row("dataset:N", title=None),
    column=alt.Column("group:N", header=alt.Header(labels=False))
).resolve_scale(
    x="independent",
)

chart = alt.vconcat(
    bp3,# bp2, bp3,
    title=alt.TitleParams(text=["Different datasets resampled and/or combined via stacked generalization.", ""], anchor="middle")
)

chart.save("chart.html")
