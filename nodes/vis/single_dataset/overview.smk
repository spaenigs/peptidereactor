import numpy as np
import pandas as pd
import altair as alt

from nodes.vis.single_dataset.scripts.utils \
    import is_struc_based

TOKEN = config["token"]

DATASET = config["dataset"]

rule overview_transform_data:
    input:
         config["metrics_dir_in"] + "{metric}.csv"
    output:
         temp(f"data/temp/{TOKEN}/{{metric}}.dp.data")
    run:
         df = pd.read_csv(input[0], index_col=0)
         df_medians = df.apply(np.median).to_frame("median")

         group = lambda enc: \
             "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:6]

         df_medians["group"] = [group(x) for x in df_medians.index]
         df_medians["metric"] = wildcards.metric
         df_medians["type"] = \
             ["structure based" if is_struc_based(e) else "sequence based" for e in df_medians["group"]]

         df_medians.to_csv(output[0])

rule collect_data:
    input:
         expand(f"data/temp/{TOKEN}/{{metric}}.dp.data",
                metric=["f1", "mcc", "precision", "recall", "sens", "spec"])
    output:
         config["html_dir_out"] + f"overview/all_metrics_dp.json"
    run:
         df_res = pd.DataFrame()
         for p in list(input):
             df_res = pd.concat([df_res, pd.read_csv(p, index_col=0)])

         df_res.to_json(output[0], orient="records")

rule create_overview_chart:
    input:
         config["html_dir_out"] + f"overview/all_metrics_dp.json"
    output:
         temp(f"data/temp/{TOKEN}/overview.json")
    run:
         url_dp = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         base = alt.Chart(url_dp).transform_aggregate(
             min_m="min(median):Q",
             max_m="max(median):Q",
             cnt_m="count(median):Q",
             groupby=["group"]
         )

         color=alt.Color("type:N", title="Type/range", scale=alt.Scale(
             domain=["sequence based", "structure based"],
             range=["#7570b3", "#d95f02"])
         )

         chart_json = alt.layer(
              base.mark_area(interpolate="monotone", opacity=0.2).encode(
                 x=alt.X('group:N', axis=alt.Axis(grid=True)),
                 y=alt.Y('min_m:Q'),
                 y2=alt.Y2('max_m:Q'),
                 color=color
             ),
             base.mark_line(interpolate="monotone", opacity=0.6).encode(
                 x=alt.X('group:N'),
                 y=alt.Y('max_m:Q'),
                 color=color
             ),
             base.mark_circle(color="red", opacity=1.0).encode(
                 x=alt.X('group:N'),
                 y=alt.Y('max_m:Q'),
                 color=color,
                 size=alt.Size("cnt_m:Q", title="Encodings per group"),
                 tooltip="cnt_m:Q"
             )
         ).facet(
             column=alt.Column("type:N", title=None),
             row=alt.Row("metric:N", title=None)
         ).resolve_scale(
             x='independent'
         ).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
