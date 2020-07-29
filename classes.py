import pandas as pd
import altair as alt

from modlamp.core import read_fasta

# dataset = "hiv_protease"
# dataset = "ace_vaxinpad"
dataset = "npp_profet"

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

base = alt.Chart(dfm)

height = 15 * len(set(source.index))

middle = base.mark_text().encode(
    y=alt.Y('range:O', axis=None),
    text=alt.Text('range:Q'),
).properties(height=height)

left = base.mark_bar().transform_filter(
    "datum.Class == 'Positive'"
).encode(
    y=alt.Y('range:O', axis=None),
    x=alt.X('Length:Q', sort=alt.Sort("descending"))
).properties(height=height)

right = base.mark_bar().transform_filter(
    "datum.Class == 'Negative'"
).encode(
    y=alt.Y('range:O', axis=None),
    x=alt.X('Length:Q')
).properties(height=height)

# TODO if too big change view

alt.concat(left, middle, right, spacing=5).save("chart.html")