import pandas as pd
import numpy as np
import altair as alt

import yaml
import re
import joblib

from nodes.vis.single_dataset.scripts.utils import cluster

TOKEN = config["token"]

with open(config["nemenyi_in"]) as f:
    nm = yaml.safe_load(f)
    CD = nm["cd"]

DOMAIN, RANGE = ["critical different", "no difference"], ["black", "gainsboro"]

rule transform_heat_map_data:
    input:
         config["crit_diff_in"]
    output:
         temp(f"data/temp/{TOKEN}/heatmap.csv")
    run:
         df_cd = pd.read_csv(input[0], index_col=0)

         row_indices = cluster(df_cd.values, axis=0)
         col_indices = cluster(df_cd.values, axis=1)
         heatmap_data = df_cd.iloc[row_indices, col_indices]

         x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
         source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "cd": heatmap_data.values.ravel()})
         source["cd_cat"] = source["cd"].apply((lambda x: "critical different" if np.abs(x) > CD else "no difference"))

         e1, e2 = [], []
         for n, s in source.iterrows():
             e1 += [df_cd.index[int(s["y"])]]
             e2 += [df_cd.columns[int(s["x"])]]

         source["Encoding1"], source["Encoding2"] = e1, e2

         source.to_csv(output[0])

rule transform_bar_chart_data:
    input:
         config["metrics_dir_in"] + "f1.csv",
         config["crit_diff_in"]
    output:
         temp(f"data/temp/{TOKEN}/barchart.csv")
    run:
         df_f1 = pd.read_csv(input[0], index_col=0)
         df_cd = pd.read_csv(input[1], index_col=0)

         dfm_count = df_f1.apply(np.mean)\
             .groupby(by=lambda x: "psekraac" if "lambda-corr" in x or "g-gap" in x else x[:6])\
             .count()\
             .to_frame("count")

         dfm_count["group"] = dfm_count.index

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

         dfm_count.to_csv(output[0])

rule transform_dots_chart_data:
    input:
         config["crit_diff_in"]
    output:
         temp(f"data/temp/{TOKEN}/dotschart.csv")
    run:
         df_cd = pd.read_csv(input[0], index_col=1)

         dataset = re.findall("data/(.*?)/", output[0])[0]

         dots_data = df_cd.values[np.triu_indices_from(df_cd.values, k=1)]
         df_dots = pd.DataFrame({"x": dots_data, "y": [dataset] * len(dots_data)})

         df_dots.to_csv(output[0])

rule create_heat_map_chart:
    input:
         f"data/temp/{TOKEN}/heatmap.csv"
    output:
         temp(f"data/temp/{TOKEN}/heatmap.joblib")
    run:
         source = pd.read_csv(input[0], index_col=0)

         hm = alt.Chart(source).mark_rect().encode(
             x=alt.X('x:O', axis=alt.Axis(title="Encoding 1", labels=False, ticks=False)),
             y=alt.X('y:O', axis=alt.Axis(title="Encoding 2", labels=False, ticks=False)),
             color=alt.Color(
                 "cd_cat:N",
                 scale=alt.Scale(domain=DOMAIN, range=RANGE),
                 legend=alt.Legend(title="CD (p<0.01)")
             ),
             tooltip=["Encoding1", "Encoding2", "cd"]
         ).properties(
             height=600, width=600
         )

         hm = alt.hconcat(
             hm,
             title=alt.TitleParams(text="All (left) vs. within group (right)", anchor="middle")
         )

         joblib.dump(hm, output[0])

rule create_bar_chart:
    input:
         f"data/temp/{TOKEN}/barchart.csv"
    output:
         temp(f"data/temp/{TOKEN}/barchart.joblib")
    run:
         dfm_count = pd.read_csv(input[0], index_col=0)

         bars1 = alt.Chart(dfm_count).mark_bar(color=RANGE[1]).encode(
             x='cd_count_max:Q',
             y="group:O",
             tooltip=["cd_count", "cd_count_max"]
         )

         bars2 = alt.Chart(dfm_count).mark_bar(color=RANGE[0]).encode(
             x=alt.X('cd_count:Q', title="# critical different"),
             y=alt.Y("group:O", title="Encoding group"),
             tooltip=["cd_count", "cd_count_max"]
         )

         bars = (bars1 + bars2).properties(height=740)

         joblib.dump(bars, output[0])

rule create_dots_chart:
    input:
         f"data/temp/{TOKEN}/dotschart.csv"
    output:
         temp(f"data/temp/{TOKEN}/dotschart.joblib")
    run:
         df_dots = pd.read_csv(input[0], index_col=0)

         dots = alt.Chart(df_dots).mark_circle(color=RANGE[0]).encode(
             x=alt.X("binned_cd:Q", title="Critical difference (binned)"),
             y=alt.Y("y:N", axis=alt.Axis(title=None)),
             size=alt.Size("count(binned_cd):Q", legend=None),
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

         joblib.dump(dots, output[0])

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

         chart = (hm & dots) | bars

         alt.data_transformers.disable_max_rows()

         with open(output[0], "w") as f:
             f.write(chart.to_json(indent=None))
             f.flush()
