import pandas as pd
import altair as alt

df = pd.read_csv("data/hiv_protease/benchmark/metrics/f1.csv", index_col=0)
df["fold"] = df.index

dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
              var_name="Encoding", value_name="F1")


selection = alt.selection_interval(encodings=["x"])
selection2 = alt.selection_interval(encodings=["x"])

first = \
    alt.Chart(dfm).mark_point().encode(
        x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
        y="average(F1)",
        opacity=alt.condition(selection, alt.value(1.0), alt.value(0.1)),
    ).add_selection(
        selection
    ).properties(width=600)

second = \
    alt.Chart(dfm).mark_point().encode(
        x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
        y="average(F1)",
    ).add_selection(
        selection2
    ).transform_filter(
        selection
    ).properties(width=600)

third = \
    alt.Chart(dfm).mark_boxplot().encode(
        x=alt.X("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
        y=alt.Y("F1"),
        opacity=alt.condition(selection2, alt.value(1.0), alt.value(0.1)),
    ).transform_filter(
        selection
    ).properties(width=600)

chart = alt.vconcat(
    # alt.Chart(dfm).transform_aggregate(
    #     mean_f1='average(F1)',
    #     groupby=["Encoding"],
    # ).mark_bar().add_selection(
    #     selection
    # ).encode(
    #     x=alt.X("mean_f1:N", bin=alt.Bin(extent=[0.0, 1.0], step=0.01)),
    #     y=alt.Y("count():Q")
    # ).properties(width=1000),
    # alt.Chart(dfm).mark_bar().add_selection(
    #     selection
    # ).encode(
    #     x=alt.X("group:Q"),
    #     y=alt.Y("count():Q")
    # )
    first, second, third
)


# import json
# import pprint
# js = json.loads(chart.to_json())
# js["datasets"] = {}
# pprint.pprint(js)
chart.save("bp2.html")