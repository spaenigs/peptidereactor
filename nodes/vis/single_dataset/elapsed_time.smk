import altair as alt
import pandas as pd

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

rule create_elapsed_time_chart:
    input:
         config["html_dir_out"] + f"elapsed_time/elapsed_time.json"
    output:
         temp(f"data/temp/{TOKEN}/elapsed_time.json")
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
                 alt.Color("category:N", title="Meta node type"),
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

         chart_json = alt.hconcat(
             c,
             title=alt.TitleParams(
                 text=["Elapsed time", ""], anchor="middle"
             ),
             config=alt.Config(
                 legend=alt.LegendConfig(titleFontSize=12, labelFontSize=12),
                 axisY=alt.AxisConfig(titleFontSize=12, titleFontWeight="normal")
             )
         ).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()