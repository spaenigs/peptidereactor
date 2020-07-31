from itertools import zip_longest
from modlamp.core import read_fasta

import Bio.Data.IUPACData as iupac_data
import altair as alt
import pandas as pd

import joblib

TOKEN = config["token"]

DATASET = config["dataset"]

rule transform_class_ratio_total_data:
    input:
         config["fasta_in"],
         config["classes_in"]
    output:
         config["html_dir_out"] + f"aac/class_ratio_total_data.json"
    run:
         seqs, names = read_fasta(input[0])
         with open(input[1]) as f:
             classes = [int(l.rstrip()) for l in f.readlines()]

         total_pos, total_neg = 0, 0
         for name, [seq, class_] in zip(names, zip(seqs, classes)):
             if class_ == 1:
                 total_pos += 1
             else:
                 total_neg += 1

         dfm_total = pd.DataFrame({"class": ["Positive", "Negative"],
                                   "count": [total_pos, total_neg]})

         dfm_total.to_json(output[0], orient="records")

rule transform_class_ratio_data:
    input:
         config["fasta_in"],
         config["classes_in"]
    output:
         config["html_dir_out"] + f"aac/class_ratio_data.json"
    run:
         seqs, names = read_fasta(input[0])
         with open(input[1]) as f:
             classes = [int(l.rstrip()) for l in f.readlines()]

         res_pos, res_neg = {}, {}
         for name, [seq, class_] in zip(names, zip(seqs, classes)):
             l = len(seq)
             if class_ == 1:
                 if l in res_pos:
                     res_pos[l] += 1
                 else:
                     res_pos[l] = 1
             else:
                 if l in res_neg:
                     res_neg[l] += 1
                 else:
                     res_neg[l] = 1

         df_pos = pd.DataFrame(res_pos.values(), index=list(res_pos.keys()), columns=["count"])\
             .sort_index()

         df_neg = pd.DataFrame(res_neg.values(), index=list(res_neg.keys()), columns=["count"])\
             .sort_index()

         source = df_pos.join(df_neg, lsuffix="_left")
         source.fillna(0, inplace=True)
         source.columns = ["Positive", "Negative"]
         source["seq_len_cat"] = source.index

         dfm = pd.melt(source, id_vars=["seq_len_cat"], value_vars=list(source.columns)[:-1],
                       var_name="class", value_name="count")

         dfm.to_json(output[0], orient="records")

rule class_ratio:
    input:
         config["html_dir_out"] + f"aac/class_ratio_data.json",
         config["html_dir_out"] + f"aac/class_ratio_total_data.json"
    output:
         temp(f"data/temp/{TOKEN}/class_ratio_chart.joblib")
    run:
         url_cr, url_total = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], ""), \
             DATASET + "/" + input[1].replace(config["html_dir_out"], "")

         filter_pos, filter_neg = \
             "datum.class == 'Positive'",  "datum.class == 'Negative'"

         dfm_total = pd.read_json(input[1])

         base_c1, base_c2 = alt.Chart(url_total), alt.Chart(url_cr)
         range_ = [0, dfm_total["count"].max()]
         height = 450

         d, r = ["Negative", "Positive"], ["#7570b3", "#d95f02"]
         color = alt.Color("class:N", scale=alt.Scale(domain=d, range=r))

         c1 = alt.hconcat(
             base_c1.mark_bar().transform_filter(filter_neg).encode(
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_),
                     sort=alt.Sort("descending")
                 ),
                 y=alt.Y('seq_len_cat:N', axis=None),
                 color=color
             ),
             base_c1.mark_bar().transform_filter(filter_pos).encode(
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_)
                 ),
                 y=alt.Y('seq_len_cat:N', axis=None),
                 color=color
             ),
             title=alt.TitleParams(text="Total count of sequences per class", anchor="middle"),
             spacing=0
         )

         c2 = alt.hconcat(
             base_c2.mark_bar().transform_filter(filter_neg).encode(
                 y=alt.Y('seq_len_cat:O', title="Sequence length"),
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_),
                     sort=alt.Sort("descending")
                 ),
                 color=color
             ),
             base_c2.mark_bar().transform_filter(filter_pos).encode(
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_)
                 ),
                 y=alt.Y('seq_len_cat:N', axis=None),
                 color=color,
             ),
             title=alt.TitleParams(
                 text="Count of sequences per length", anchor="middle"
             ),
             spacing=0
         )

         joblib.dump(alt.vconcat(c1, c2, spacing=8), output[0])

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

         df = pd.DataFrame({"Negative": neg, "Positive": pos})
         df["Negative"] = df["Negative"].apply(lambda c: (c / df["Negative"].sum()) * 100)
         df["Positive"] = df["Positive"].apply(lambda c: (c / df["Positive"].sum()) * 100)
         df["aa"] = df.index

         df_melted = pd.melt(df, id_vars=["aa"], value_vars=list(df.columns)[:-1],
                             var_name="class", value_name="count")

         df_melted.to_json(output[0], orient="records")

rule aac:
    input:
         config["html_dir_out"] + f"aac/aac_data.json"
    output:
         temp(f"data/temp/{TOKEN}/aac_chart.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         chart = alt.Chart(url).mark_bar().encode(
             x=alt.X("class:N", axis=None),
             y=alt.Y("count:Q", title="Relative count (%)", scale=alt.Scale(domain=[0.0, 100.0])),
             color=alt.Color(
                 "class:N",
                 scale=alt.Scale(
                     domain=["Negative", "Positive"],
                     range=["#7570b3", "#d95f02"])),
             column=alt.Column("aa:N", title="Amino acid"),
             tooltip="count:Q"
         ).properties(
             width=20
         )

         joblib.dump(chart, output[0])

rule collect_aac:
    input:
         f"data/temp/{TOKEN}/class_ratio_chart.joblib",
         f"data/temp/{TOKEN}/aac_chart.joblib"
    output:
         temp(f"data/temp/{TOKEN}/amino_acid_comp.json")
    run:
         cr, aac = joblib.load(input[0]), joblib.load(input[1])

         chart_json = (cr & aac).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
