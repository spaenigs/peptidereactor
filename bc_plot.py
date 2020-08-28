import pandas as pd
import altair as alt

import numpy as np
from plotly.validators.parcoords import _labelangle

dataset = "hiv_protease"
# dataset = "ace_vaxinpad"


def order(x):
    print(x)
    if "ensemble" in x and "resampled" in x:
        return f"z2_{x}"
    elif "ensemble" in x:
        return f"z1_{x}"
    elif "resampled" in x:
        return f"{x}"
    else:
        return f"{x}"

# df = pd.read_csv(f"bc_data_{dataset}.csv")
df = pd.read_csv(f"/home/spaenigs/Downloads/bc_data_{dataset}.csv")
df["encoding"] = df["encoding"].apply(lambda x: order(x))
print(df.head())

import re
from scipy import stats

# for encoding_pair in df["group"].unique():
#     e1, e2 = re.findall("'(.*?)', '(.*?)'", encoding_pair)[0]
#     e1_resampled_f1s = df.loc[df["encoding"].str.contains("resampled"), :].loc[df["encoding"].str.contains(e1), "f1"].values
#     e1_non_resampled_f1s = df.loc[~df["encoding"].str.contains("resampled"), :].loc[df["encoding"].str.contains(e1), "f1"].values
#     _, e1_resampled_vs_non_pv = stats.ttest_rel(e1_non_resampled_f1s, e1_resampled_f1s)
#     print(e1_resampled_vs_non_pv)

bp = alt.Chart(df).mark_boxplot().encode(
    x=alt.X("encoding:N", title=None),
    y='f1:Q',
    color="resampled:N",
    column=alt.Column("group:N", header=alt.Header(labels=False))
).resolve_scale(
    x="independent"
)

scatter = alt.Chart(df).mark_point(filled=True, size=60).encode(
    x=alt.X("encoding:N", title=None, sort="y"),
    y='median(f1):Q',
    color="resampled:N",
    column=alt.Column("group:N", header=alt.Header(labels=False)),
    tooltip=["encoding:N", 'median(f1):Q']
).resolve_scale(
    x="independent"
)

(bp & scatter).save("chart.html")
