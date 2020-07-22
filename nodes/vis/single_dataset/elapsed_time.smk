import altair as alt
import pandas as pd

TOKEN = config["token"]

rule create_elapsed_time_chart:
    input:
         config["benchmark_csv_in"]
    output:
         temp(f"data/temp/{TOKEN}/elapsed_time.json")
    run:
         df_time = pd.read_csv(input[0], index_col=0)
         df_time["meta_rule"] = df_time.index

         d, rc, rs = \
             ["encoding", "benchmark", "utils"], \
             ["#1b9e77", "#d95f02", "#7570b3"], \
             ["cross", "square", "triangle-right"]

         c = alt.Chart(df_time).mark_point(opacity=1.0, size=50).encode(
             x=alt.X("meta_rule:N", sort="-y", title="Meta rule"),
             y=alt.Y("time_in_hours:Q", scale=alt.Scale(type='log'),  title="Hours (log scale)"),
             color=alt.Color("category:N"),
             shape=alt.Shape("category:N"),
             tooltip=["meta_rule", "time_in_hours:Q", "category:N"]
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