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

         dfm_total = pd.DataFrame({"class": ["1", "0"],
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
         source.columns = ["1", "0"]
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
             "datum.class == '1'",  "datum.class == '0'"

         dfm_total = pd.read_json(input[1])

         base_c1, base_c2 = alt.Chart(url_total), alt.Chart(url_cr)
         range_ = [0, dfm_total["count"].max()]
         height = 450

         d, r = ["0", "1"], ["#7b3294", "#008837"]
         color = alt.Color("class:N", title="Class label", scale=alt.Scale(domain=d, range=r))

         c1 = alt.hconcat(
             base_c1.mark_bar().transform_filter(filter_neg).encode(
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_),
                     sort=alt.Sort("descending")
                 ),
                 y=alt.Y('seq_len_cat:N', axis=None),
                 color=color,
                 tooltip="count:Q"
             ),
             base_c1.mark_bar().transform_filter(filter_pos).encode(
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_)
                 ),
                 y=alt.Y('seq_len_cat:N', axis=None),
                 color=color,
                 tooltip="count:Q"
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
                 color=color,
                 tooltip="count:Q"
             ),
             base_c2.mark_bar().transform_filter(filter_pos).encode(
                 x=alt.X('count:Q',
                     title="Count",
                     scale=alt.Scale(domain=range_)
                 ),
                 y=alt.Y('seq_len_cat:N', axis=None),
                 color=color,
                 tooltip="count:Q"
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

         df = pd.DataFrame({"0": neg, "1": pos})
         df["0"] = df["0"].apply(lambda c: (c / df["0"].sum()) * 100)
         df["1"] = df["1"].apply(lambda c: (c / df["1"].sum()) * 100)
         df["aa"] = df.index

         df_melted = pd.melt(df, id_vars=["aa"], value_vars=list(df.columns)[:-1],
                             var_name="class", value_name="count")

         df_melted["anno"] = 0

         df_melted.to_json(output[0], orient="records")

rule aac:
    input:
         config["html_dir_out"] + f"aac/aac_data.json"
    output:
         temp(f"data/temp/{TOKEN}/aac_chart.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         def anno_layer(val):
             return alt.Chart().mark_rule(
                 opacity=0.5,
                 strokeWidth=0.3
             ).encode(
                 y="anno_tmp:Q"
             ).transform_calculate(
                 anno_tmp=alt.datum.anno + val
             )

         chart = alt.layer(
             alt.Chart().mark_bar().encode(
                 x=alt.X("class:N", axis=None),
                 y=alt.Y("count:Q", title="Relative count (%)", scale=alt.Scale(domain=[0.0, 100.0])),
                 color=alt.Color(
                     "class:N",
                     scale=alt.Scale(
                         domain=["0", "1"],
                         range=["#7b3294", "#008837"]
                     )
                 ),
                 tooltip="tt:N"
             ).transform_calculate(
                 tt="join(['~', round(datum.count), '%'], '')"
             ).properties(
                 width=18
             ),
             anno_layer(20),
             anno_layer(40),
             anno_layer(60),
             anno_layer(80)
         ).facet(
             column=alt.Column("aa:N", title="Amino acid"),
             data=url
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

         chart_json = alt.vconcat(
             cr, aac,
             center=True,
             title=alt.TitleParams(
                 text=[
                     "Overall class distribution (top), sequence length distribution (middle) and",
                     "amino acid ratios relative to the respective class (bottom).",
                     "",
                     ""
                 ],
                 anchor="middle"
             ),
             config=alt.Config(
                 view=alt.ViewConfig(width=600),
                 legend=alt.LegendConfig(titleFontSize=12, labelFontSize=12),
                 axis=alt.AxisConfig(titleFontSize=12, titleFontWeight="normal")
             )
         ).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
