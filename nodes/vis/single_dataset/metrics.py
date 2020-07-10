import pandas as pd
import altair as alt


def metrics_chart(df_f1, df_mcc):

    def chart(df_metrics, metric):

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

        def melt_and_annotate(df, metric):
            df["fold"] = df.index
            dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
                          var_name="Encoding", value_name="Value")
            dfm["metric"] = metric
            return dfm

        dfm = melt_and_annotate(df_metrics, metric)
        dfm = dfm.drop("fold", axis=1).groupby(by=["Encoding"]).median().reset_index()
        dfm.sort_values(by="Value", ascending=False, inplace=True)
        dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["Encoding"]]

        nr_top_encodings = 20

        dfc = df_metrics[dfm["Encoding"][:nr_top_encodings]].copy()
        dfm1 = melt_and_annotate(dfc, metric)
        dfm1["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm1["Encoding"]]

        hline_data = pd.DataFrame({'a': [dfm["Value"][nr_top_encodings]]})
        anno_struc = dfm.loc[dfm["type"] == "structure based", :]

        d, r = ["sequence based", "structure based"], ["#7570b3", "#d95f02"]

        scatter = alt.layer(
            alt.Chart(dfm).mark_circle().encode(
                x=alt.X(
                    "Encoding:N", sort="-y",
                    axis=alt.Axis(labels=False, ticks=False, title=None)
                ),
                y=alt.Y(
                    "Value:Q", title=metric,
                    axis=alt.Axis(grid=False), scale=alt.Scale(domain=[0.0, 1.0])
                ),
                size=alt.value(20),
                opacity=alt.condition(
                    alt.datum.Value >= dfm["Value"][nr_top_encodings],
                    alt.value(1.0),
                    alt.value(0.3)
                ),
                color=alt.Color("type:N", scale=alt.Scale(domain=d, range=r)),
                tooltip="Encoding:N"
            ).properties(
                width=600
            ),
            alt.Chart(hline_data).mark_rule(color="grey", strokeDash=[1, 1], opacity=0.5).encode(y="a:Q"),
            alt.Chart(anno_struc).mark_rule(opacity=0.3).encode(
                x=alt.X("Encoding:N", sort="-y"),
                y="Value:Q",
                color="type:N"
            )
        )

        bp = alt.Chart(dfm1).mark_boxplot().encode(
            x=alt.X("Encoding", title=None),
            y=alt.Y("Value", title=None, scale=alt.Scale(domain=[0.0, 1.0])),
            color=alt.Color("type:N", scale=alt.Scale(domain=d, range=r)),
        ).properties(
            width=600
        )

        return scatter | bp

    c1 = chart(df_f1, "F1")
    c2 = chart(df_mcc, "MCC")

    return c1 & c2


# def get_overview_data(df, metric):
#     df["fold"] = df.index
#     dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
#                   var_name="Encoding", value_name="Value")
#
#     def is_struc_based(e):
#         if "asa" in e:
#             return True
#         elif "delaunay" in e:
#             return True
#         elif "disorder" in e:
#             return True
#         elif "hull" in e:
#             return True
#         elif "qsar" in e:
#             return True
#         elif "sse" in e:
#             return True
#         elif "ta" in e:
#             return True
#         else:
#             return False
#
#     dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["Encoding"]]
#     dfm["metric"] = metric
#     return dfm

# test2 = test2.add_selection(
#     selection4
# ).properties(width=600, height=800)

# selection = alt.selection_interval(encodings=["x"])
# selection2 = alt.selection_interval(encodings=["x"])
# selection3 = alt.selection_multi(fields=['type'], bind='legend')
# selection4 = alt.selection_single()
#
# first = \
#     alt.Chart(dfm).mark_point().encode(
#         x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
#         y="average(Value)",
#         # opacity=alt.condition(selection, alt.value(1.0), alt.value(0.1)),
#         color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
#         shape=alt.condition(select_dropdown, alt.Shape("metric:N"), alt.Color("white")),
#         tooltip="Encoding:N"
#     ).add_selection(
#         selection
#     ).add_selection(
#         selection3
#     ).transform_filter(
#         selection3
#     ).properties(width=600)
#
# second = \
#     alt.Chart(dfm).mark_point().encode(
#         x=alt.X("Encoding", sort="-y"), #, axis=alt.Axis(labels=False, ticks=False, title=None)),
#         y="average(Value)",
#         color="type:N",
#         shape="metric:N",
#     ).add_selection(
#         selection2
#     ).transform_filter(
#         selection
#     ).properties(width=600)
#
# third = \
#     alt.Chart(dfm).mark_boxplot().encode(
#         y=alt.Y("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
#         x=alt.X("Value"),
#         opacity=alt.condition(selection2, alt.value(1.0), alt.value(0.1)),
#         color="type:N",
#         shape="metric:N",
#     ).transform_filter(
#         selection
#     ).add_selection(selection4).properties(width=600, height=800)

# # chart = alt.hconcat(alt.vconcat(
# #     # alt.Chart(dfm).transform_aggregate(
# #     #     mean_f1='average(F1)',
# #     #     groupby=["Encoding"],
# #     # ).mark_bar().add_selection(
# #     #     selection
# #     # ).encode(
# #     #     x=alt.X("mean_f1:N", bin=alt.Bin(extent=[0.0, 1.0], step=0.01)),
# #     #     y=alt.Y("count():Q")
# #     # ).properties(width=1000),
# #     # alt.Chart(dfm).mark_bar().add_selection(
# #     #     selection
# #     # ).encode(
# #     #     x=alt.X("group:Q"),
# #     #     y=alt.Y("count():Q")
# #     # )
# #     first, second,
# #         alt.Chart(dfm).mark_bar().encode(
# #             x=alt.X("F1", bin=alt.Bin(maxbins=20)),
# #             y="count()",
# #             opacity=alt.condition(selection4, alt.value(1.0), alt.value(0.1)),
# #         ).transform_filter(selection4)
# # ), third
# # )
#
# df = pd.read_csv("data/hiv_protease/benchmark/similarity/seq_vs_str/diversity.csv", index_col=0)
#


# chart = alt.vconcat(
#     alt.hconcat(alt.vconcat(first, second), third),
#     hm
# ).properties(
#     title="HIV protease"
# ).configure_title(
#     anchor="middle",
#     fontSize=26,
#     color="purple"
# )

# import json
# import pprint
# js = json.loads(chart.to_json())
# js["datasets"] = {}
# pprint.pprint(js)
# chart.save("bp2.html")


# encodings = sorted(list(set(dfm["Encoding"])))
#
# dropdown = alt.binding_select(options=[None] + encodings)
# dd_selector = alt.selection_single(fields=["Encoding"], bind=dropdown, name="Value of")
#
# slider = alt.binding_range(min=-1.0, max=1.0, step=0.05, name='cutoff:')
# selector = alt.selection_single(fields=['avg_value'], bind=slider, init={'avg_value': 0.0})
#
# return alt.Chart(dfm).transform_joinaggregate(
#     groupby=["Encoding", "metric"],
#     avg_value="average(Value)"
# ).mark_point().encode(
#     x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=True, ticks=True)),
#     y=alt.Y("average(Value)", scale=alt.Scale(domain=[-1.0, 1.0])),
#     color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
#     shape="metric:N",
#     tooltip="Encoding:N"
# ).add_selection(
#     selector, dd_selector
# ).transform_filter(
#     alt.datum.avg_value > selector.avg_value
# ).transform_filter(
#     dd_selector
# ).properties(
#     width=1430,
#     height=300
# )

# input_dropdown = alt.binding_select(options=["F1", "MCC", None])
# select_dropdown = alt.selection_single(fields=["metric"], bind=input_dropdown, name="Value of",
#                                        init={"metric": "F1"})
# selection = alt.selection_interval(encodings=["x"], init={"x": [dfm["Encoding"][0], dfm["Encoding"][20]]})

# scatter = alt.Chart(dfm).transform_joinaggregate(
#     groupby=["Encoding", "metric"],
#     avg_value="average(Value)"
# ).mark_point().encode(
#     x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
#     y=alt.Y("avg_value:O", scale=alt.Scale(domain=[-1.0, 1.0])),
#     color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
#     shape="metric:N",
#     tooltip="Encoding:N"
# ).transform_filter(
#     alt.datum.avg_value > 0.8
# ).transform_filter(
#     select_dropdown
# ).add_selection(
#     select_dropdown
# ).add_selection(
#     selection
# ).properties(
#     width=600,
#     height=300
# )

# bp = alt.Chart(dfm).mark_boxplot().encode(
#     y=alt.Y("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
#     x=alt.X("Value"),
#     # opacity=alt.condition(selection2, alt.value(1.0), alt.value(0.1)),
#     color="type:N",
#     shape="metric:N"
# ).transform_filter(
#     select_dropdown
# ).transform_filter(
#     selection
# ).properties(
#     width=600,
#     height=300
# )

# return alt.hconcat(scatter, bp)
