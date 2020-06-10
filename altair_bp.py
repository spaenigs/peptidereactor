import pandas as pd
import altair as alt

df = pd.read_csv("data/hiv_protease/benchmark/metrics/f1.csv", index_col=0)
df["fold"] = df.index

dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
              var_name="Encoding", value_name="F1")


def is_struc_based(e):
    if "asa" in e:
        return True
    elif "delaunay" in e:
        return True
    elif "disorder" in e:
        return True
    elif "hull" in e:
        return True
    elif "qsar" in e:
        return True
    elif "sse" in e:
        return True
    elif "ta" in e:
        return True
    else:
        return False


dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["Encoding"]]

selection = alt.selection_interval(encodings=["x"])
selection2 = alt.selection_interval(encodings=["x"])
selection3 = alt.selection_multi(fields=['type'], bind='legend')
selection4 = alt.selection_single()

first = \
    alt.Chart(dfm).mark_point().encode(
        x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
        y="average(F1)",
        opacity=alt.condition(selection, alt.value(1.0), alt.value(0.1)),
        color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
        tooltip="Encoding:N"
    ).add_selection(
        selection
    ).add_selection(
        selection3
    ).transform_filter(
        selection3
    ).properties(width=600)

second = \
    alt.Chart(dfm).mark_point().encode(
        x=alt.X("Encoding", sort="-y"), #, axis=alt.Axis(labels=False, ticks=False, title=None)),
        y="average(F1)",
        color="type:N"
    ).add_selection(
        selection2
    ).transform_filter(
        selection
    ).properties(width=600)

third = \
    alt.Chart(dfm).mark_boxplot().encode(
        y=alt.Y("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
        x=alt.X("F1"),
        opacity=alt.condition(selection2, alt.value(1.0), alt.value(0.1)),
        color="type:N"
    ).transform_filter(
        selection
    ).add_selection(selection4).properties(width=600, height=800)

chart = alt.hconcat(alt.vconcat(
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
    first, second,
        alt.Chart(dfm).mark_bar().encode(
            x=alt.X("F1", bin=alt.Bin(maxbins=20)),
            y="count()",
            opacity=alt.condition(selection4, alt.value(1.0), alt.value(0.1)),
        ).transform_filter(selection4)
), third
)


# import json
# import pprint
# js = json.loads(chart.to_json())
# js["datasets"] = {}
# pprint.pprint(js)
chart.save("bp2.html")