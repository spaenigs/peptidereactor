import altair as alt
import pandas as pd
import numpy as np

import joblib

from nodes.vis.single_dataset.scripts.utils \
    import is_struc_based

TOKEN = config["token"]

DATASET = config["dataset"]

rule dump_elapsed_time_data:
    input:
         config["benchmark_csv_in"]
    output:
         config["html_dir_out"] + f"elapsed_time/elapsed_time.json"
    run:
         df_time = pd.read_csv(input[0], index_col=0)
         df_time["meta_rule"] = df_time.index
         df_time.to_json(output[0], orient="records")

rule dump_time_vs_performance_data:
    input:
        config["metrics_dir_in"] + "f1.csv",
        config["benchmark_csv_in"]
    output:
        config["html_dir_out"] + f"elapsed_time/time_vs_performance.json"
    run:
        df = pd.read_csv(input[0], index_col=0)
        df_time = pd.read_csv(input[1], index_col=0)

        df_medians = df.apply(np.median).to_frame("median")

        group = lambda enc: \
            "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:6]

        df_medians["group"] = [group(x) for x in df_medians.index]
        df_medians = df_medians.groupby("group").median()
        df_medians["group"] = df_medians.index
        df_medians.index = range(0, df_medians.shape[0])
        df_medians["type"] = \
            ["structure based" if is_struc_based(e) else "sequence based" for e in df_medians["group"]]

        df_time["meta_rule"] = df_time.index
        df_time["group"] = [group(x) for x in df_time["meta_rule"]]

        dur_tert_search = \
            df_time.loc[df_time["meta_rule"]=="tertiary_structure_search", "s"][0]/60/60

        times = []
        for n, s in df_medians.iterrows():
            sg = s["group"]
            df_tmp = df_time.loc[[g in sg and sg.startswith(g) for g in df_time["group"]], :]
            if sg == "aainde":
                df_tmp = df_time.loc[df_time["meta_rule"] == "encoding_aaindex", :]
            elif sg == "dist_f":
                df_tmp = df_time.loc[df_time["meta_rule"] == "distance_frequency", :]
            if not df_tmp.empty:
                times += [df_tmp["s"][0]]
            else:
                print(f"Could find {sg} in elapsed time df. Skipping.")
                df_medians.drop(df_medians.loc[df_medians["group"] == sg, :].index, inplace=True)

        df_medians["time"] = [t/60/60 for t in times]
        df_medians["time"] = df_medians\
            .apply(lambda s: s["time"] + dur_tert_search if is_struc_based(s.group) else s["time"], axis=1)\
            .to_list()

        df_medians.to_json(output[0], orient="records")


rule create_time_vs_performance_chart:
    input:
         config["html_dir_out"] + f"elapsed_time/time_vs_performance.json"
    output:
         temp(f"data/temp/{TOKEN}/time_vs_performance.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         c = alt.Chart(url).mark_point(filled=True).encode(
             x=alt.X("time:Q", title="Hours (log scale)", scale=alt.Scale(type="log")),
             y=alt.Y("median:Q", title="Median performance of encoding groups"),
             color=alt.Color(
                 "type:N",
                 title="Encoding type",
                 scale=alt.Scale(
                     domain=["sequence based", "structure based"],
                     range=["#7570b3", "#d95f02"]
                 )
             ),
             tooltip=["group:N", "time:Q"]
         )

         text = c.mark_text(dy=10, size=9).encode(
             text="group:N"
         )

         joblib.dump((c + text).interactive(), output[0])

rule create_elapsed_time_chart:
    input:
         config["html_dir_out"] + f"elapsed_time/elapsed_time.json"
    output:
         temp(f"data/temp/{TOKEN}/elapsed_time.joblib")
    run:
         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         d, rc = \
             ["encoding", "benchmark", "utils"], \
             ["#1b9e77", "#d95f02", "#7570b3"]

         selection = alt.selection_single(on="mouseover", fields=["category"])

         c = alt.Chart(url).mark_point(filled=True, size=50).encode(
             x=alt.X(
                 "meta_rule:N", sort="-y", title=None, axis=alt.Axis(grid=True)
             ),
             y=alt.Y(
                 "time_in_hours:Q",
                 title="Hours (log scale)",
                 scale=alt.Scale(type='log'),
                 axis=alt.Axis(grid=False)
             ),
             color=alt.condition(
                 selection,
                 alt.Color("category:N", title="Meta node type", scale=alt.Scale(scheme="category10", reverse=True)),
                 alt.value("grey")
             ),
             opacity=alt.condition(selection, alt.value(1.0), alt.value(0.0)),
             tooltip=["meta_rule:N", "time_in_hours:Q", "category:N"]
         ).transform_calculate(
             time_in_hours="datum.s/60/60"
         ).add_selection(
             selection
         ).properties(
             width=1200
         )

         joblib.dump(c, output[0])

rule concat_time_charts:
    input:
         f"data/temp/{TOKEN}/time_vs_performance.joblib",
         f"data/temp/{TOKEN}/elapsed_time.joblib"
    output:
         temp(f"data/temp/{TOKEN}/elapsed_time.json")
    run:
         time_vs_perf = joblib.load(input[0])
         elapsed_time = joblib.load(input[1])

         chart_json = alt.vconcat(
             time_vs_perf, elapsed_time,
             title=alt.TitleParams(
                 text=[
                     "Median of group performance vs. elapsed time (encodings only) (top) and",
                     "elapsed time for all meta nodes (bottom).",
                     ""
                 ]
                 , anchor="middle"
             ),
             center=True,
             config=alt.Config(
                 legend=alt.LegendConfig(titleFontSize=12, labelFontSize=12),
                 axisX=alt.AxisConfig(titleFontSize=12, titleFontWeight="normal"),
                 axisY=alt.AxisConfig(titleFontSize=12, titleFontWeight="normal")
             )
         ).resolve_scale(
             color="independent"
         ).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()