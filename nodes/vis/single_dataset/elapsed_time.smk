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

         d, rc, rs = \
             ["encoding", "benchmark", "utils"], \
             ["#1b9e77", "#d95f02", "#7570b3"], \
             ["cross", "square", "triangle-right"]

         c = alt.Chart(url).mark_point(opacity=1.0, size=50).encode(
             x=alt.X("meta_rule:N", sort="-y", title="Meta rule"),
             y=alt.Y("time_in_hours:Q", scale=alt.Scale(type='log'),  title="Hours (log scale)"),
             color=alt.Color("category:N"),
             shape=alt.Shape("category:N"),
             tooltip=["meta_rule:N", "time_in_hours:Q", "category:N"]
         ).transform_calculate(
             time_in_hours="datum.s/60/60"
         ).add_selection(
             alt.selection_single()
         ).properties(
             width=900
         )

         chart_json = alt.hconcat(c, title=alt.TitleParams(
             text=["Elapsed time", ""], anchor="middle"
         )).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()