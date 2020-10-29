import pandas as pd
import altair as alt

from glob import glob
from modlamp.core import read_fasta

import re

from more_itertools import intersperse

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
print(df_res.head())

df_res["desc"] = "a" * 15 + "\n" + "b" * 15 + "\n" + "a" * 15 + "\n" + "b" * 15 + "\n" + "a" * 15 + "\n" + "b" * 15 + "\n"
df_res["ref"] = "www.heiderlab.de"


def wrap_text(text):
    word = text.strip()
    words = word.split(" ")
    if len(words) > 1:
        wrapped_sentence = ""
        tmp_sentence = ""
        for w in words:
            if len(tmp_sentence + w) <= 15:
                tmp_sentence += w + " "
            else:
                wrapped_sentence += tmp_sentence + "\n"
                tmp_sentence = w + " "
        return wrapped_sentence + tmp_sentence
    else:
        return "".join(intersperse("\n", word, n=15))


desc_col, ref_col = [], []
for ds in df_res.dataset:
    try:
        with open(f"data/{ds}/README.md") as f:
            col1, col2 = [], []
            for h1, h2 in re.findall("\|(.*?)\|(.*?)\|", f.read()):
                col1 += [h1]
                col2 += [h2]
            df_res.loc[df_res.dataset == ds, "desc"] = wrap_text(col1[-1])
            df_res.loc[df_res.dataset == ds, "ref"] = wrap_text(col2[-1])
            # desc_col += [col1[-1]]
            # ref_col += [col2[-1]]
    except Exception:
        pass

print(df_res.head())

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

ds = alt.Chart().mark_text(
    fontSize=12,
    lineBreak="\n",
).encode(
    text="dataset:N"
).properties(
    height=80,
    width=125
)

ss = alt.Chart().mark_text(
    fontSize=12
).encode(
    text="seq_size:O",
).properties(
    height=80,
    width=125
)

t = alt.Chart().mark_text(
    dy=-30,
    fontSize=12,
    lineBreak="\n",
).encode(
    text="desc:N"
).properties(
    height=80,
    width=125
)

ref = alt.Chart().transform_calculate(
    url="" + alt.datum.ref
).mark_text(
    dy=-30,
    fontSize=12,
    lineBreak="\n"
).encode(
    text="ref:N",
    href="url:N",
    tooltip="url:N"
).properties(
    height=80,
    width=125
)

table = alt.hconcat(ds, ss, t, ref, data=df_res).transform_filter(
   alt.datum.dataset == "hiv_ddi"
)

df_tsne = pd.read_json("data/multiple_datasets/vis/md_tsne/tsne_data.json")

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

alt.vconcat(alt.hconcat(scatter, tsnec), table).save("chart.html", vegalite_version="4.17.0")
# table.save("chart.html", vegalite_version="4.17.0")

