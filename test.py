import pandas as pd
import altair as alt

from glob import glob
from modlamp.core import read_fasta

import re

in_ = glob("data/*/misc/benchmark/benchmark.csv")

res = []
for p in list(in_):
    dataset = re.findall("data/(.*?)/", p)[0]
    df = pd.read_csv(p, index_col=0)
    res += [[dataset, df["s"].sum() / 60 / 60]]

df_time = pd.DataFrame(res, columns=["dataset", "hours"])

fin_ = glob("data/*/seqs.fasta")

res = []
for p in list(fin_):
    dataset = re.findall("data/(.*?)/", p)[0]
    seqs, _ = read_fasta(p)
    res += [[dataset, len(seqs)]]

df_len = pd.DataFrame(res, columns=["dataset", "seq_size"])

df_res = pd.merge(df_time, df_len, on="dataset")
# df_res["bio_field"] = df_res.dataset.apply(lambda ds: ds[:3])

selection = alt.selection_single(
    fields=["dataset"],
    init={"dataset": "hiv_ddi"},
    empty="none"
)

scatter = alt.Chart(
    df_res,
    title="Multiple datasets"
).mark_point(filled=True, size=60).encode(
    x=alt.X(
        "hours:Q",
        title="log(Computation time) (h)",
        scale=alt.Scale(type='log'),
        axis=alt.Axis(grid=False, titleFontWeight="normal")
    ),
    y=alt.Y(
        "seq_size:Q",
        title="log(# of sequences)",
        sort="-x",
        scale=alt.Scale(type='log'),
        axis=alt.Axis(grid=False, titleFontWeight="normal")
    ),
    tooltip="dataset:N",
    color=alt.condition(selection, alt.value("#4C78A8"), alt.value("lightgrey"))
).add_selection(
    selection
).properties(
    height=250,
    width=250
).interactive()

df_tsne = pd.read_json("data/multiple_datasets/vis/md_tsne/tsne_data.json")
# df_tsne["bio_field"] = df_tsne.dataset.apply(lambda ds: ds[:3])

# print(df_tsne.head())

x_min, x_max, y_min, y_max = -100, 100, -100, 100

scatterc = alt.Chart().mark_circle(
    size=10,
    color="#fdc086"
).encode(
    x=alt.X(
        "x:Q",
        title="tSNE-1",
        axis=alt.Axis(grid=False, titleFontWeight="normal"),
        scale=alt.Scale(domain=[x_min, x_max])
    ),
    y=alt.Y(
        "y:Q",
        title="tSNE-2",
        axis=alt.Axis(grid=False, titleFontWeight="normal"),
        scale=alt.Scale(domain=[y_min, y_max])
    )
).transform_filter(
    selection
)

hullc = alt.Chart().mark_line(
    color="#386cb0",
    strokeDash=[5, 3],
    strokeWidth=1
).encode(
    x="x:Q",
    y="y:Q",
    order="order:O"
).transform_filter(
    alt.datum.hull_vertex == True
).transform_filter(
    selection
)

textc = alt.Chart().mark_text().encode(
    x="x:Q",
    y="y:Q",
    text="text:N"
).transform_calculate(
    text="join(['area=', round(datum.area)], '')"
).transform_filter(
    (alt.datum.x == 0) and (alt.datum.y == 0)
).transform_filter(
    selection
)

tsnec = alt.layer(
    scatterc, hullc, textc,
    data=df_tsne,
    title="Single dataset"
).properties(height=250, width=250)

alt.hconcat(scatter, tsnec).save("chart.html")

