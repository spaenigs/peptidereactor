import pandas as pd
import numpy as np
import altair as alt

import yaml
import re
import joblib

from nodes.vis.single_dataset.scripts.utils import cluster, is_struc_based

TOKEN = config["token"]

DATASET = config["dataset"]

CRIT_DIFF_DIR_IN = config["crit_diff_dir_in"]

with open(CRIT_DIFF_DIR_IN + "nemenyi.yaml") as f:
    nm = yaml.safe_load(f)
    CD = nm["cd"]

DOMAIN, RANGE = ["critical (p<0.01)", "non-critical"], ["black", "gainsboro"]

rule transform_heat_map_data:
    input:
         CRIT_DIFF_DIR_IN + "diff_matrix.csv"
    output:
         config["html_dir_out"] + f"crit_diff/heatmap.json"
    run:
         df_cd = pd.read_csv(input[0], index_col=0)

         row_indices = cluster(df_cd.values, axis=0)
         col_indices = cluster(df_cd.values, axis=1)
         heatmap_data = df_cd.iloc[row_indices, col_indices]

         x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
         source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "cd": heatmap_data.values.ravel()})
         source["cd_cat"] = source["cd"].apply((lambda x: "critical (p<0.01)" if np.abs(x) > CD else "non-critical"))

         source["Encoding1"] = source["x"].apply(lambda i: heatmap_data.columns[i])
         source["Encoding2"] = source["y"].apply(lambda i: heatmap_data.index[i])

         source.to_json(output[0], orient="records")

rule transform_bar_chart_data:
    input:
         config["metrics_dir_in"] + "f1.csv",
         CRIT_DIFF_DIR_IN + "diff_matrix.csv"
    output:
         config["html_dir_out"] + f"crit_diff/barchart.json"
    run:
         df_f1 = pd.read_csv(input[0], index_col=0)
         df_cd = pd.read_csv(input[1], index_col=0)

         dfm_count = df_f1.apply(np.mean)\
             .groupby(by=lambda x: "psekraac" if "lambda-corr" in x or "g-gap" in x else x[:6])\
             .count()\
             .to_frame("count")

         dfm_count["group"] = dfm_count.index

         dfm_count["type"] = \
             ["structure based" if is_struc_based(e) else "sequence based" for e in dfm_count["group"]]

         for i in dfm_count.index:

             if i == "psekraac":
                 indices_filtered = [ii for ii in df_cd.index if "g-gap" in ii or "lambda" in ii]
             else:
                 indices_filtered = [ii for ii in df_cd.index if i == ii[:6]]

             df_cd_sub = df_cd.loc[indices_filtered, indices_filtered]
             df_le = df_cd_sub.le(-CD)
             df_ge = df_cd_sub.ge(CD)

             r = sum(df_le.values[np.triu_indices_from(df_le.values, k=1)]) + \
                 sum(df_ge.values[np.triu_indices_from(df_ge.values, k=1)])

             dfm_count.loc[i, "cd_count"] = r
             dfm_count.loc[i, "cd_count_max"] = \
                 len(indices_filtered) * (len(indices_filtered) - 1) / 2

         dfm_count = dfm_count.loc[dfm_count["count"] > 1, :]
         dfm_count.to_json(output[0], orient="records")

rule transform_dots_chart_data:
    input:
         CRIT_DIFF_DIR_IN + "diff_matrix.csv"
    output:
         config["html_dir_out"] + f"crit_diff/dotschart.json"
    run:
         df_cd = pd.read_csv(input[0], index_col=1)

         dataset = re.findall("data/(.*?)/", output[0])[0]

         dots_data = df_cd.values[np.triu_indices_from(df_cd.values, k=1)]
         df_dots = pd.DataFrame({"x": dots_data, "y": [dataset] * len(dots_data)})

         df_dots.to_json(output[0], orient="records")

rule create_heat_map_chart:
    input:
         config["html_dir_out"] + f"crit_diff/heatmap.json"
    output:
         temp(f"data/temp/{TOKEN}/heatmap.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         hm = alt.Chart(
             url,
             title="All vs. all encodings"
         ).mark_rect().encode(
             x=alt.X('x:O', axis=alt.Axis(title="Encoding 1", labels=False, ticks=False)),
             y=alt.X('y:O', axis=alt.Axis(title="Encoding 2", labels=False, ticks=False)),
             color=alt.Color(
                 "cd_cat:N",
                 scale=alt.Scale(domain=DOMAIN, range=RANGE),
                 legend=alt.Legend(title="Difference")
             ),
             tooltip=["Encoding1:N", "Encoding2:N", "cd:Q"]
         ).properties(
             height=600, width=600
         )

         joblib.dump(hm, output[0])

rule create_bar_chart:
    input:
         config["html_dir_out"] + f"crit_diff/barchart.json"
    output:
         temp(f"data/temp/{TOKEN}/barchart.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         bars1 = alt.Chart(url).mark_bar(color=RANGE[1]).encode(
             x='cd_count_max:Q',
             y="group:O",
             tooltip=["cd_count:Q", "cd_count_max:Q"]
         )

         bars2 = alt.Chart(url).mark_bar(color=RANGE[0]).encode(
             x=alt.X('cd_count:Q', title="Count"),
             y=alt.Y("group:O", title=None),
             tooltip=["cd_count:Q", "cd_count_max:Q"]
         )

         bars = alt.layer(
             bars1, bars2,
             title=alt.TitleParams(
                 text="Within encoding groups (group size > 1)",
                 anchor="middle"
             )
         ).facet(
             row=alt.Row(
                 "type:N",
                 title=None,
                 header=alt.Header(labels=False)
             )
         ).resolve_scale(
             y="independent"
         )

         joblib.dump(bars, output[0])

rule create_dots_chart:
    input:
         config["html_dir_out"] + f"crit_diff/dotschart.json"
    output:
         temp(f"data/temp/{TOKEN}/dotschart.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         dots = alt.Chart(url).mark_circle(color=RANGE[0]).encode(
             x=alt.X("binned_cd:Q", title="Difference"),
             y=alt.Y("y:N", axis=alt.Axis(title=None)),
             size=alt.Size(
                 "count(binned_cd):Q",
                 title="Bin size",
                 # legend=alt.LegendConfig(orient="bottom")
             ),
             color=alt.condition(np.abs(alt.datum.binned_cd) >= CD,
                                 alt.value(RANGE[0]),
                                 alt.value(RANGE[1])),
             tooltip=["count(binned_cd):Q", "binned_cd:Q"]
         ).transform_bin(
             "binned_cd", "x", bin=alt.Bin(extent=[-300, 300], step=30)
         ).properties(
             width=600,
             height=100
         )

         df_lines = pd.DataFrame({
             "x": [-CD, CD],
             "cat": ["negative", "positive"]
         })

         lines = alt.Chart(
             df_lines,
             title="Binned difference for all encoding pairings"
         ).mark_rule(strokeDash=[5, 2]).encode(
             x="x:Q",
             color=alt.Color(
                 "cat:N",
                 title="Critical difference",
                 scale=alt.Scale(
                     domain=["negative", "positive"],
                     range=["#d8b365", "#5ab4ac"]
                 ),
                 # legend=alt.LegendConfig(orient="bottom")
             )
         )

         joblib.dump(lines + dots, output[0])

rule concat_cd_charts:
    input:
         f"data/temp/{TOKEN}/heatmap.joblib",
         f"data/temp/{TOKEN}/barchart.joblib",
         f"data/temp/{TOKEN}/dotschart.joblib"
    output:
         temp(f"data/temp/{TOKEN}/critical_difference.json")
    run:
         hm = joblib.load(input[0])
         bars = joblib.load(input[1])
         dots = joblib.load(input[2])

         chart = alt.vconcat(
             (hm | bars), dots,
             title=alt.TitleParams(
                 text=[
                     "Statistical comparison (critical difference) of model performance per fold,",
                     "for all vs. all encodings, within encoding group and",
                     "the binned counts distribution of differences for all comparisions."
                     "",
                     "",
                     ""
                 ],
                 anchor="middle"
             ),
             config=alt.Config(
                 legend=alt.LegendConfig(titleFontSize=12, labelFontSize=12),
                 axis=alt.AxisConfig(titleFontSize=12, titleFontWeight="normal")
             )
         ).resolve_scale(size="independent", color="independent")

         alt.data_transformers.disable_max_rows()

         with open(output[0], "w") as f:
             f.write(chart.to_json(indent=None))
             f.flush()
