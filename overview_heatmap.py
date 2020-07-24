from functools import reduce
from glob import glob

import pandas as pd
import altair as alt
import numpy as np

import re


def compute_median_f1(path):
    dataset = re.findall("data/(.*?)/", path)[0]
    return pd\
        .read_csv(path, index_col=0)\
        .apply(np.median)\
        .to_frame(dataset)


paths = glob("data/*/benchmark/metrics/f1.csv")

df_res = pd.concat(
    axis=1,
    objs=reduce(
        lambda dfs, p: dfs + [compute_median_f1(p)],
        paths, [])).\
    sort_index().\
    fillna(0.0).\
    T

x, y = np.meshgrid(df_res.columns, df_res.index)

source = pd.DataFrame({"Encoding": x.ravel(),
                       "Dataset": y.ravel(),
                       "F1": df_res.values.ravel()
                       })

l = np.hstack((np.array([df_res.index]).T, df_res.values))
l = sorted(l, key=lambda x: [*x[1:]], reverse=True)
sorted_datasets = np.array(l)[:, 0]

alt.Chart(source).mark_rect().encode(
    x='Encoding:O',
    y=alt.Y(
        'Dataset:N',
        sort=alt.Sort(alt.SortArray(sorted_datasets))
    ),
    color='F1:Q',
    tooltip=["Encoding", "Dataset", "F1"]
).properties(
    width=1700
).save("chart1.html")