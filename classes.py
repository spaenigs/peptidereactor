import pandas as pd
import altair as alt

from vega_datasets import data
from modlamp.core import read_fasta

# dataset = "hiv_protease"
dataset = "ace_vaxinpad"
# dataset = "npp_profet"
# dataset = "afp_amppred"

seqs, names = read_fasta(f"data/{dataset}/seqs.fasta")
with open(f"data/{dataset}/classes.txt") as f:
    classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

seq_tuples_pos = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)) if tup[1] == 1)
seq_tuples_neg = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)) if tup[1] == 0)

pos = [len(tup[0]) for tup in seq_tuples_pos.values()]
neg = [len(tup[0]) for tup in seq_tuples_neg.values()]

df1 = pd.DataFrame({"seq_lengths": pos, "class": ["positive" for i in range(0, len(pos))]})\
    .groupby(by="seq_lengths")\
    .count()

df2 = pd.DataFrame({"seq_lengths": neg, "class": ["negative" for i in range(0, len(neg))]})\
    .groupby(by="seq_lengths")\
    .count()

source = df1.join(df2, lsuffix="_left")
source.fillna(0, inplace=True)
source.columns = ["Positive", "Negative"]
source["range"] = source.index

dfm = pd.melt(source, id_vars=["range"], value_vars=list(source.columns)[:-1],
              var_name="Class", value_name="Length")

dfm.groupby("Class").apply(lambda df: (df["Length"]/df["Length"].sum())*100).to_frame("ratio")

base = alt.Chart(dfm)

range = [0, dfm["Length"].sort_values().unique()[-1]]

height = 450

left = base.mark_bar().transform_filter(
    "datum.Class == 'Positive'"
).encode(
    x=alt.X('Length:Q', scale=alt.Scale(domain=range), sort=alt.Sort("descending")),
    y=alt.Y('range:O'),
    color="Class:N"
).properties(height=height)

right = base.mark_bar().transform_filter(
    "datum.Class == 'Negative'"
).encode(
    y=alt.Y('range:O', axis=None),
    x=alt.X('Length:Q', scale=alt.Scale(domain=range)),
    color="Class:N"
).properties(height=height)

cnts = dfm.groupby(by="Class")["Length"].sum().to_frame("nr")
cnts["Class"] = cnts.index

range2 = [0, cnts["nr"].max()]
print(range2)

bar_left = alt.Chart(cnts).mark_bar().transform_filter(
    "datum.Class == 'Positive'"
).encode(
    y=alt.Y('Class:N'),
    x=alt.X('nr:Q', scale=alt.Scale(domain=range2), sort=alt.Sort("descending")),
    color="Class:N"
).properties(height=20)

bar_right = alt.Chart(cnts).mark_bar().transform_filter(
    "datum.Class == 'Negative'"
).encode(
    y=alt.Y('Class:N', axis=None),
    x=alt.X('nr:Q', scale=alt.Scale(domain=range2)),
    color="Class:N"
).properties(height=20)

c1 = alt.hconcat(left, right, title=alt.TitleParams(text="All (left) vs. within group (right)", anchor="middle"))
c2 = alt.hconcat(bar_left, bar_right, title=alt.TitleParams(text="All (left) vs. within group (right)", anchor="middle"))

alt.vconcat(c1, c2, spacing=0).save("chart.html")


