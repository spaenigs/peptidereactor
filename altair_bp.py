import pandas as pd
import altair as alt
import seaborn as sns
import numpy as np

two_charts_template = """
<!DOCTYPE html>
<html>
<head>
  <script src="https://cdn.jsdelivr.net/npm/vega@{vega_version}"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@{vegalite_version}"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@{vegaembed_version}"></script>
</head>
<body>

<center><p><b>{dataset}</b></p></center> 
 
<p><u>General overview</u></p>
<div id="vis1"></div>

<p><u>Metrics</u></p>
<div id="vis2"></div>

<p><u>Similarity</u></p>
<div id="vis3"></div>

<script type="text/javascript">
  vegaEmbed('#vis1', {spec1}).catch(console.error);
  vegaEmbed('#vis2', {spec2}).catch(console.error);
  vegaEmbed('#vis3', {spec3}).catch(console.error);
</script>
</body>
</html>
"""


def main():

    alt.data_transformers.disable_max_rows()

    df_f1 = pd.read_csv("data/hiv_protease/benchmark/metrics/f1.csv", index_col=0)
    df_mcc = pd.read_csv("data/hiv_protease/benchmark/metrics/mcc.csv", index_col=0)
    dfm = pd.concat([get_overview_data(df_f1, "F1"), get_overview_data(df_mcc, "MCC")])

    df_div = pd.read_csv("data/hiv_protease/benchmark/similarity/seq_vs_str/diversity.csv", index_col=0)
    df_phi = pd.read_csv("data/hiv_protease/benchmark/similarity/seq_vs_str/phi.csv", index_col=0)

    with open('bp2.html', 'w') as f:
        f.write(two_charts_template.format(
            vega_version=alt.VEGA_VERSION,
            vegalite_version=alt.VEGALITE_VERSION,
            vegaembed_version=alt.VEGAEMBED_VERSION,
            dataset="HIV protease",
            spec1=get_overview_chart(dfm).to_json(indent=None),
            spec2=get_metrics_chart(dfm).to_json(indent=None),
            spec3=get_similarity_chart(df_div, df_phi).to_json(indent=None)
        ))


def get_similarity_chart(df_div, df_phi):

    def get_clustering(df):
        heatmap_data = sns.clustermap(df.values).data2d
        x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
        source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "Similarity": heatmap_data.values.ravel()})
        e1, e2 = [], []
        for n, s in source.iterrows():
            idx = heatmap_data.index[int(s["y"])]
            col = heatmap_data.columns[int(s["x"])]
            e1 += [df.index[idx]]
            e2 += [df.columns[col]]
        source["Encoding1"], source["Encoding2"] = e1, e2
        return source

    def get_heatmap_chart(source, title):
        hm = alt.Chart(source).mark_rect().encode(
            x=alt.X('x:O', title="Encoding 2", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.X('y:O', title="Encoding 1", axis=alt.Axis(labels=False, ticks=False)),
            color='Similarity:Q',
            tooltip=["Encoding1", "Encoding2", "Similarity"]
        ).add_selection(
            alt.selection_single()
        ).properties(
            title=title,
            width=710,
            height=600
        )
        return hm

    return alt.hconcat(
        get_heatmap_chart(get_clustering(df_div), "Diversity"),
        get_heatmap_chart(get_clustering(df_phi), "Phi correlation")
    )


def get_metrics_chart(dfm):
    input_dropdown = alt.binding_select(options=["F1", "MCC", None])
    select_dropdown = alt.selection_single(fields=["metric"], bind=input_dropdown, name="Value of",
                                           init={"metric": "F1"})
    selection = alt.selection_interval(encodings=["x"])

    scatter = alt.Chart(dfm).mark_point().encode(
        x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
        y=alt.Y("average(Value)", scale=alt.Scale(domain=[-1.0, 1.0])),
        color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
        shape="metric:N",
        tooltip="Encoding:N"
    ).transform_filter(
        select_dropdown
    ).add_selection(
        select_dropdown
    ).add_selection(
        selection
    ).properties(
        width=600,
        height=300
    )

    bp = alt.Chart(dfm).mark_boxplot().encode(
        y=alt.Y("Encoding", sort=alt.EncodingSortField(field="F1", op='mean')),
        x=alt.X("Value"),
        # opacity=alt.condition(selection2, alt.value(1.0), alt.value(0.1)),
        color="type:N",
        shape="metric:N"
    ).transform_filter(
        select_dropdown
    ).transform_filter(
        selection
    ).properties(
        width=600,
        height=300
    )

    return alt.hconcat(scatter, bp)


def get_overview_chart(dfm):
    slider = alt.binding_range(min=-1.0, max=1.0, step=0.05, name='cutoff:')
    selector = alt.selection_single(fields=['avg_value'], bind=slider, init={'avg_value': 0.0})
    return alt.Chart(dfm).transform_joinaggregate(
        groupby=["Encoding", "metric"],
        avg_value="average(Value)"
    ).mark_point().encode(
        x=alt.X("Encoding", sort="-y", axis=alt.Axis(labels=True, ticks=True)),
        y=alt.Y("average(Value)", scale=alt.Scale(domain=[-1.0, 1.0])),
        color=alt.Color("type:N", scale=alt.Scale(domain=["structure based", "sequence based"], range=["#e7ba52", "#1f77b4"])),
        shape="metric:N",
        tooltip="Encoding:N"
    ).add_selection(
        selector
    ).transform_filter(
        alt.datum.avg_value > selector.avg_value
    ).properties(
        width=1430,
        height=300
    )


def get_overview_data(df, metric):
    df["fold"] = df.index
    dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
                  var_name="Encoding", value_name="Value")
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
    dfm["metric"] = metric
    return dfm


if __name__ == '__main__':
    main()






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