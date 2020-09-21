from glob import glob
from functools import reduce

import pandas as pd
import numpy as np
import altair as alt

import re

from nodes.vis.single_dataset.scripts.utils \
    import cluster, is_struc_based


def is_imbalanced(path):
    with open(path) as f:
        classes = [int(i.rstrip()) for i in f.readlines()]
        return sum(classes) / len(classes)


def compute_median_f1(path):
    dataset = re.findall("data/(.*?)/", path)[0]
    return pd\
        .read_csv(path, index_col=0)\
        .apply(np.median)\
        .to_frame(dataset)


inp = glob("data/*/benchmark/metrics/")

data_out = "data.json"
out = "chart.html"

paths = [p + "f1.csv" for p in inp]

df_res = pd.DataFrame()
df_group_size = pd.DataFrame()
for p in paths:

    df = pd.read_csv(p, index_col=0)
    df_medians = df.apply(np.median).to_frame("median")
    group = lambda enc: \
        "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:6]
    df_medians["group"] = [group(x) for x in df_medians.index]

    df_tmp = df_medians.groupby(by="group").max()
    df_tmp.columns = [re.findall("data/(.*?)/", p)[0]]
    df_res = pd.concat([df_res, df_tmp], axis=1)

    df_group_size_tmp = df_medians.groupby(by="group").count()
    df_group_size_tmp["encoding"] = df_group_size_tmp.index
    df_group_size_tmp.columns = ["size", "encoding"]
    df_group_size_tmp["type"] = \
        ["structure based" if is_struc_based(e) else "sequence based" for e in df_group_size_tmp["encoding"]]
    df_group_size = pd.concat([df_group_size, df_group_size_tmp], axis=0)

x, y = np.meshgrid(df_res.columns, df_res.index)

source = pd.DataFrame({"Dataset": x.ravel(),
                       "Encoding": y.ravel(),
                       "F1": df_res.values.ravel()
                       })

source["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in source["Encoding"]]

source["is_imbalanced"] = source["Dataset"].apply(lambda ds: is_imbalanced(f"data/{ds}/classes.txt"))
source.sort_values(by="is_imbalanced", inplace=True)
names_imbalanced = source["Dataset"].unique()

fields = sorted(source["Dataset"].apply(lambda x: x[:3]).unique())
fields_dict = dict((k, v) for k, v in zip(fields, range(1, len(fields)+1)))
source["bio_field"] = source["Dataset"].apply(lambda ds: ds[:3])

sorted_encodings = source\
    .groupby(by="Encoding")\
    .mean().sort_values("F1", ascending=False)\
    .index.to_list()

hcharts = []
for b in sorted(source["bio_field"].unique()):
    vcharts = []
    for t in source["type"].unique():
        df_tmp = source.loc[source["bio_field"] == b, :].loc[source["type"] == t, :]
        df_tmp.sort_values(by="is_imbalanced", inplace=True)
        vcharts += [alt.Chart(df_tmp).mark_rect(strokeWidth=5.0).encode(
                y=alt.Y(
                    'Encoding:N',
                    title=None if not b.startswith("ace") else \
                        "Sequence-based encodings" if t == "sequence based" else "Structure-based encodings",
                    axis=alt.Axis(labels=False if not b.startswith("ace") else True,
                                  ticks=False if not b.startswith("ace") else True),
                    sort=alt.Sort(alt.SortArray(sorted(sorted_encodings)))
                ),
                x=alt.Y(
                    'Dataset:N',
                    title=None,
                    axis=alt.Axis(labels=False if t == "sequence based" else True,
                                  ticks=False if t == "sequence based" else True),
                    sort=alt.Sort(alt.SortArray(source["Dataset"].unique()))
                ),
                color=alt.Color("F1:Q", title="F1"),
                tooltip=["Encoding:N", "Dataset:N", "F1:Q", "is_imbalanced:Q"]
            ).properties(
                height=500 if t == "sequence based" else 150
            )
        ]
    hcharts += [alt.vconcat(*vcharts, spacing=2)]

chart2 = alt.hconcat(*hcharts, spacing=5)

###

fields_dict = {"sequence based": 1, "structure based": 2}
source["type_field"] = source["type"].apply(lambda t: fields_dict[t])

vcharts = []
for t in source["type"].unique():
    hcharts = []
    for thresholds in [[0.0, 0.35], [0.35, 1.0]]:
        df_tmp = source.loc[source["type"] == t, :].loc[source["is_imbalanced"].between(*thresholds), :]
        hcharts += [alt.Chart(df_tmp).mark_rect().encode(
            y=alt.Y(
                'Encoding:N',
                title=None,
                # title=None if thresholds[1] == 1.0 else \
                #     "Sequence-based encodings" if t == "sequence based" else "Structure-based encodings",
                axis=alt.Axis(labels=False, ticks=False)
                # axis=alt.Axis(labels=False if thresholds[1] == 1.0 else True,
                #               ticks=False if thresholds[1] == 1.0 else True),
            ),
            x=alt.X(
                'Dataset:N',
                sort=alt.Sort(alt.SortArray(names_imbalanced)),
                title=None,
                # title=None if t == "sequence based" else ["Imbalanced", "datasets"] if thresholds[1] == 0.3 else [
                #     "Balanced", "datasets"],
                axis=alt.Axis(labels=False if t == "sequence based" else True,
                              ticks=False if t == "sequence based" else True)
            ),
            color=alt.Color("F1:Q", title="F1"),
            tooltip=["Encoding:N", "Dataset:N", "F1:Q", "is_imbalanced:Q"]
        ).properties(
            height=500 if t == "sequence based" else 150
        )]
    vcharts += [alt.hconcat(*hcharts, spacing=2)]

chart3 = alt.vconcat(*vcharts, spacing=2)

vcharts = []
for t in source["type"].unique():
    df_tmp = df_group_size.loc[df_group_size["type"] == t, :].groupby("encoding").mean()
    df_tmp = df_tmp.loc[[str(i).endswith(".0") for i in df_tmp["size"]], :]
    df_tmp = df_group_size.loc[df_group_size["encoding"].isin(df_tmp.index)].drop_duplicates()
    vcharts += [
        alt.Chart(df_group_size.loc[df_group_size["type"] == t, :].drop_duplicates()).mark_boxplot().encode(
            y=alt.Y(
                "encoding:N",
                title=None,
                axis=alt.Axis(labels=False, ticks=False)
            ),
            x=alt.X(
                "size:Q",
                title=None,
                axis=alt.Axis(labels=False if t == "sequence based" else True,
                              ticks=False if t == "sequence based" else True)
            )
        ).properties(
            height=500 if t == "sequence based" else 150,
            width=200
        ) + alt.Chart(df_tmp).mark_point(shape="square", filled=True).encode(
            x="size:Q",
            y="encoding:N",
            tooltip=["encoding:N", "size:Q"]
        )
    ]

chart4 = alt.vconcat(*vcharts, spacing=2).resolve_scale(
    x="shared"
)

alt.hconcat(chart2, chart4, chart3).save(out)
