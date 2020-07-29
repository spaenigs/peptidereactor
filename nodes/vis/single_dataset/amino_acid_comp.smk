from itertools import zip_longest
from modlamp.core import read_fasta

import Bio.Data.IUPACData as iupac_data
import altair as alt
import pandas as pd

TOKEN = config["token"]

DATASET = config["dataset"]

rule aac_transform_data:
    input:
         config["fasta_in"],
         config["classes_in"]
    output:
         config["html_dir_out"] + f"aac/aac_data.json"
    run:
         seqs, names = read_fasta(input[0])
         with open(input[1]) as f:
             classes = [int(l.rstrip()) for l in f.readlines()]

         neg = dict(zip_longest(iupac_data.protein_letters, [0] * 20))
         pos = dict(zip_longest(iupac_data.protein_letters, [0] * 20))

         aas = lambda c: [seq for seq, class_ in zip(seqs, classes) if class_ == c]

         for aa in "".join(aas(0)):
             neg[aa] += 1

         for aa in "".join(aas(1)):
             pos[aa] += 1

         df = pd.DataFrame({"neg": neg, "pos": pos})
         df["neg"] = df["neg"].apply(lambda c: (c / df["neg"].sum()) * 100)
         df["pos"] = df["pos"].apply(lambda c: (c / df["pos"].sum()) * 100)
         df["aa"] = df.index

         df_melted = pd.melt(df, id_vars=["aa"], value_vars=list(df.columns)[:-1],
                             var_name="class", value_name="count")

         df_melted.to_json(output[0], orient="records")

rule aac:
    input:
         config["html_dir_out"] + f"aac/aac_data.json"
    output:
         temp(f"data/temp/{TOKEN}/amino_acid_comp.json")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         chart_json = alt.Chart(url).mark_bar().encode(
             x=alt.X("class:N", axis=None),
             y=alt.Y("count:Q", title="Relative count (%)", scale=alt.Scale(domain=[0.0, 100.0])),
             color=alt.Color(
                 "class:N",
                 scale=alt.Scale(
                     domain=["neg", "pos"],
                     range=["#7570b3", "#d95f02"])),
             column=alt.Column("aa:N", title="Amino acid"),
             tooltip="count:Q"
         ).properties(
             width=20
         ).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
