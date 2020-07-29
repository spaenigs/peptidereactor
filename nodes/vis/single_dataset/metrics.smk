import pandas as pd
import altair as alt

import joblib

from nodes.vis.single_dataset.scripts.utils \
    import is_struc_based

TOKEN = config["token"]

DATASET = config["dataset"]
NR_TOP_ENCODINGS = 20

rule metrics_transform_data:
    input:
         config["metrics_dir_in"] + "{metric}.csv"
    output:
         config["html_dir_out"] + f"metrics/{{metric}}_scatter_plot_data.json",
         config["html_dir_out"] + f"metrics/{{metric}}_box_plot_data.json"
    run:
         def melt_and_annotate(df, metric):
             df["fold"] = df.index
             dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
                           var_name="Encoding", value_name="Value")
             dfm["metric"] = metric
             return dfm

         df_metrics = pd.read_csv(input[0], index_col=0)

         scatter_plot_data = melt_and_annotate(df_metrics, wildcards.metric)
         scatter_plot_data = scatter_plot_data.drop("fold", axis=1)\
             .groupby(by=["Encoding"])\
             .median()\
             .reset_index()
         scatter_plot_data.sort_values(by="Value", ascending=False, inplace=True)
         scatter_plot_data["type"] = \
             ["structure based" if is_struc_based(e) else "sequence based" for e in scatter_plot_data["Encoding"]]

         dfc = df_metrics[scatter_plot_data["Encoding"][:NR_TOP_ENCODINGS]].copy()
         box_plot_data = melt_and_annotate(dfc, wildcards.metric)
         box_plot_data["type"] = \
             ["structure based" if is_struc_based(e) else "sequence based" for e in box_plot_data["Encoding"]]

         scatter_plot_data.to_json(output[0], orient="records")
         box_plot_data.to_json(output[1], orient="records")

rule create_metrics_chart:
    input:
         config["html_dir_out"] + f"metrics/{{metric}}_scatter_plot_data.json",
         config["html_dir_out"] + f"metrics/{{metric}}_box_plot_data.json"
    output:
         config["html_dir_out"] + f"metrics/{{metric}}_hline_plot_data.json",
         config["html_dir_out"] + f"metrics/{{metric}}_anno_plot_data.json",
         config["html_dir_out"] + f"metrics/{{metric}}_vline_plot_data.json",
         temp(f"data/temp/{TOKEN}/{{metric}}.joblib")
    run:
         scatter_plot_data = pd.read_json(input[0])
         box_plot_data = pd.read_json(input[1])

         hline_data = pd.DataFrame({'a': [scatter_plot_data["Value"][NR_TOP_ENCODINGS]]})
         anno_struc = scatter_plot_data.loc[scatter_plot_data["type"] == "structure based", :]
         vline_data = pd.DataFrame({'y': [0.0]})

         metric = wildcards.metric

         # adopt url, such that it can be found from server
         repl_path = \
             lambda path: DATASET + "/" + path.replace(config["html_dir_out"], "")

         url_scatter, url_box, url_hline, url_anno, url_vline = \
             input[0], input[1], output[0], output[1], output[2]

         d, r = ["sequence based", "structure based"], ["#7570b3", "#d95f02"]

         scatter = alt.layer(
             alt.Chart(repl_path(url_scatter)).mark_circle().encode(
                 x=alt.X(
                     "Encoding:N", sort="-y",
                     axis=alt.Axis(labels=False, ticks=False, title=None)
                 ),
                 y=alt.Y(
                     "Value:Q", title=metric,
                     axis=alt.Axis(grid=False), scale=alt.Scale(
                         domain=[-1.0 if metric == "mcc" else 0.0, 1.0]
                     )
                 ),
                 size=alt.value(20),
                 opacity=alt.condition(
                     alt.datum.Value >= scatter_plot_data["Value"][NR_TOP_ENCODINGS],
                     alt.value(1.0),
                     alt.value(0.3)
                 ),
                 color=alt.Color("type:N", scale=alt.Scale(domain=d, range=r)),
                 tooltip="Encoding:N"
             ).properties(
                width=600
             ),
             alt.Chart(repl_path(url_hline)).mark_rule(
                 color="grey",
                 strokeDash=[1, 1],
                 opacity=0.5).encode(
                 y="a:Q"
             ),
             alt.Chart(repl_path(url_anno)).mark_rule(opacity=0.3).encode(
                 x=alt.X("Encoding:N", sort="-y"),
                 y="Value:Q",
                 color="type:N"
             )
         )

         if metric == "mcc":
             scatter = scatter + alt.Chart(repl_path(url_vline)).mark_rule(
                 color="black",
                 opacity=0.5
             ).encode(
                 y="y:Q"
             )

         bp = alt.Chart(repl_path(url_box)).mark_boxplot().encode(
             x=alt.X("Encoding:N", title=None),
             y=alt.Y("Value:Q", title=None, scale=alt.Scale(domain=[0.0, 1.0])),
             color=alt.Color("type:N", scale=alt.Scale(domain=d, range=r)),
         ).properties(
             width=600
         )

         hline_data.to_json(url_hline, orient="records")
         anno_struc.to_json(url_anno, orient="records")
         vline_data.to_json(url_vline, orient="records")

         joblib.dump(scatter | bp, output[3])

rule concat_charts:
    input:
         expand(f"data/temp/{TOKEN}/{{metric}}.joblib",
                metric=["f1", "mcc", "precision", "recall", "sens", "spec"])
    output:
         temp(f"data/temp/{TOKEN}/metrics.json")
    run:
         charts = [joblib.load(p) for p in list(input)]

         chart_json = alt.vconcat(*charts).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
