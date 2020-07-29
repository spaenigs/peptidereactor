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
         temp(f"data/temp/{TOKEN}/{{metric}}.data")
    run:
         df = pd.read_csv(input[0], index_col=0)

         dfm_count = df.apply(np.median)\
             .groupby(by=lambda x: "psekraac" if "lambda-corr" in x or "g-gap" in x else x[:6])\
             .count()\
             .to_frame("count")
         dfm_max = df.apply(np.median)\
             .groupby(by=lambda x: "psekraac" if "lambda-corr" in x or "g-gap" in x else x[:6])\
             .max()\
             .to_frame("max_metric")

         dfm = pd.concat([dfm_max, dfm_count], axis=1)
         dfm["group"] = dfm.index
         dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["group"]]
         dfm["metric"] = wildcards.metric

         dfm.to_csv(output[0])

rule collect_data:
    input:
         expand(f"data/temp/{TOKEN}/{{metric}}.data",
                metric=["f1", "mcc", "precision", "recall", "sens", "spec"])
    output:
         config["html_dir_out"] + f"overview/all_metrics.json"
    run:
         df_res = pd.DataFrame()
         for p in list(input):
             df_res = pd.concat([df_res, pd.read_csv(p, index_col=0)])

         df_res.to_json(output[0], orient="records")

rule create_overview_chart:
    input:
         config["html_dir_out"] + f"overview/all_metrics.json"
    output:
         temp(f"data/temp/{TOKEN}/overview.json")
    run:
         url = DATASET + "/" + input[0].replace(config["html_dir_out"], "")
         base = alt.Chart(url)

         chart_json =  alt.layer(
             base.mark_circle(filled=True, size=50, opacity=0.7).encode(
                 x=alt.X("group:N", axis=alt.Axis(title=None)),
                 y=alt.Y("max_metric:Q", axis=alt.Axis(title=None), scale=alt.Scale(domain=[-0.2, 1.0])),
                 size=alt.Size("count:Q"),
                 color=alt.Color(
                     "type:N",
                     scale=alt.Scale(
                         domain=["sequence based", "structure based"],
                         range=["#7570b3", "#d95f02"])
                 ),
                 tooltip="count:Q"
             ),
             base.mark_rule(opacity=1.0).encode(
                 x=alt.X("group:N"), y="max_metric:Q", color="type:N"
             )
         ).facet(
             column=alt.Column("type:N", title=None),
             row=alt.Row("metric:N", title=None)
         ).resolve_scale(
             x='independent'
         ).to_json(
             indent=None
         )

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
