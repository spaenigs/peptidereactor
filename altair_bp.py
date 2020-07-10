import pandas as pd
import altair as alt
import seaborn as sns
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score, auc
import matplotlib.pyplot as plt

two_charts_template = """
<!DOCTYPE html>
<html>
<head>
  <script src="https://cdn.jsdelivr.net/npm/vega@{vega_version}"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@{vegalite_version}"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@{vegaembed_version}"></script>
</head>
<body>

<center>
<img src="/home/spaenigs/Downloads/peptide-reactor/peptide-reactor/peptide-reactor.005.png" alt="logo" width="7%" height="7%"/>
<hr>

<p><b>{dataset}</b></p>

</center> 
 
<p><u>General overview</u></p>
<center><div id="vis1"></div></center>

<p><u>Metrics</u></p>
<center><div id="vis2"></div></center>

<p><u>ROC</u></p>
<center><div id="vis3"></div></center>

<p><u>Similarity</u></p>
<center><div id="vis4"></div></center>

<script type="text/javascript">
  vegaEmbed('#vis1', {spec1}).catch(console.error);
  vegaEmbed('#vis2', {spec2}).catch(console.error);
  vegaEmbed('#vis3', {spec3}).catch(console.error);
  vegaEmbed('#vis4', {spec4}).catch(console.error);
</script>
</body>
</html>
"""


def main():

    alt.data_transformers.disable_max_rows()

    dataset = "hiv_protease"

    df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    df_mcc = pd.read_csv(f"data/{dataset}/benchmark/metrics/mcc.csv", index_col=0)
    dfm = pd.concat([get_overview_data(df_f1, "F1"), get_overview_data(df_mcc, "MCC")])

    df_div = pd.read_csv(f"data/{dataset}/benchmark/similarity/seq_vs_str/diversity.csv", index_col=0)
    df_phi = pd.read_csv(f"data/{dataset}/benchmark/similarity/seq_vs_str/phi.csv", index_col=0)
    df_div_1 = pd.read_csv(f"data/{dataset}/benchmark/similarity/all_vs_all/diversity.csv", index_col=0)
    df_phi_1 = pd.read_csv(f"data/{dataset}/benchmark/similarity/all_vs_all/phi.csv", index_col=0)
    source = pd.concat([get_clustering(df_div, "Diversity", "Seq vs. struc."),
                        get_clustering(df_phi, "Phi correlation", "Seq vs. struc."),
                        # get_clustering(df_div_1, "Diversity", "All vs. all"),
                        # get_clustering(df_phi_1, "Phi correlation", "All vs. all")
                        ])

    df_test = pd.read_csv("data/hiv_protease/benchmark/single/y_true_cv_ngram_a3_300.csv", index_col=0)
    df_prob = pd.read_csv("data/hiv_protease/benchmark/single/y_prob_cv_ngram_a3_300.csv", index_col=0)
    df_test1 = pd.read_csv("data/hiv_protease/benchmark/single/y_true_cv_delaunay_total_distance.csv", index_col=0)
    df_prob1 = pd.read_csv("data/hiv_protease/benchmark/single/y_prob_cv_delaunay_total_distance.csv", index_col=0)
    res = pd.concat(
        [get_roc_data(df_test, df_prob, "ngram_a3_300"), get_roc_data(df_test1, df_prob1, "delaunay_total_distance")])

    with open('bp2.html', 'w') as f:
        f.write(two_charts_template.format(
            vega_version=alt.VEGA_VERSION,
            vegalite_version=alt.VEGALITE_VERSION,
            vegaembed_version=alt.VEGAEMBED_VERSION,
            dataset=dataset,
            spec1=get_overview_chart(df_f1, df_mcc).to_json(indent=None),
            spec2=get_metrics_chart(df_f1, df_mcc).to_json(indent=None),
            spec3=get_roc_chart(res).to_json(indent=None),
            spec4=get_similarity_chart(source).to_json(indent=None)
        ))


def get_roc_data(df_test, df_prob, encoding):
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    for i in range(df_test.shape[0]):
        y_true = df_test.iloc[i, :].dropna().values
        y_pred = df_prob.iloc[i, :].dropna().values
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(roc_auc_score(y_true, y_pred))
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    df = pd.DataFrame({"mean_fpr": mean_fpr, "mean_tpr": mean_tpr, "tprs_lower": tprs_lower, "tprs_upper": tprs_upper})
    df["Encoding"] = encoding
    return df


def get_roc_chart(res):

    return alt.Chart(res).mark_line().encode(
        x="mean_fpr",
        y="mean_tpr",
        color="Encoding:N"
    )


def get_clustering(df, measure, comparision):
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
    source["Similarity measure"] = measure
    source["Comparision"] = comparision
    return source


def get_similarity_chart(source):

    selection_highlight = alt.selection_single(on="click", resolve="global") #, clear="mousemove", resolve="global")

    hm = alt.Chart(source).mark_rect().encode(
        x=alt.X('x:O', title="Encoding 2", axis=alt.Axis(labels=False, ticks=False)),
        y=alt.X('y:O', title="Encoding 1", axis=alt.Axis(labels=False, ticks=False)),
        color="Similarity:Q",
        opacity=alt.condition(selection_highlight, alt.value(1.0), alt.value(0.6)),
        tooltip=["Encoding1", "Encoding2", "Similarity"]
    ).add_selection(
        selection_highlight
    ).properties(
        width=720,
        height=600
    ).facet(
        column="Similarity measure:N",
        row="Comparision:N"
    )

    return hm


def get_metrics_chart(df_f1, df_mcc):

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

    def get_overview_data(df, metric):
        df = df.copy()
        df["fold"] = df.index
        dfm = pd.melt(df, id_vars=["fold"], value_vars=list(df.columns)[:-1],
                      var_name="Encoding", value_name="Value")
        dfm["metric"] = metric
        return dfm

    def get_chart(df, metric):

        dfm = get_overview_data(df, metric)
        dfm = dfm.drop("fold", axis=1).groupby(by=["Encoding"]).median().reset_index()
        # dfm["Encoding"] = dfm.index
        dfm.sort_values(by="Value", ascending=False, inplace=True)
        dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["Encoding"]]

        nr_top_encodings = 20

        hline_data = pd.DataFrame({'a': [dfm["Value"][nr_top_encodings]]})
        anno_struc = dfm.loc[dfm["type"] == "structure based", :]

        scatter = alt.layer(
            alt.Chart(dfm).mark_circle().encode(
                x=alt.X("Encoding:N", sort="-y", axis=alt.Axis(labels=False, ticks=False, title=None)),
                y=alt.Y("Value:Q", title=metric, axis=alt.Axis(grid=False), scale=alt.Scale(domain=[0.0, 1.0])),
                size=alt.value(20),
                opacity=alt.condition(alt.datum.Value >= dfm["Value"][nr_top_encodings], alt.value(1.0), alt.value(0.3)),
                color=alt.Color("type:N", scale=alt.Scale(domain=["sequence based", "structure based"],
                                                          range=["#7570b3", "#d95f02"])),
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

        dfm = get_overview_data(df[dfm["Encoding"][:nr_top_encodings]], metric)
        dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["Encoding"]]

        bp = alt.Chart(dfm).mark_boxplot().encode(
            x=alt.X("Encoding", title=None),
            y=alt.Y("Value", title=None, scale=alt.Scale(domain=[0.0, 1.0])),
            # color=alt.value("red"),
            color=alt.Color("type:N", scale=alt.Scale(domain=["sequence based", "structure based"],
                                                      range=["#7570b3", "#d95f02"])),
        ).properties(
            width=600
        )

        return scatter | bp

    c1 = get_chart(df_f1, "F1")
    c2 = get_chart(df_mcc, "MCC")

    return c1 & c2


def get_overview_chart(df_f1, df_mcc):

    def is_struc_based(e):
        if "asa" in e:
            return True
        elif "delaun" in e:
            return True
        elif "disor" in e:
            return True
        elif "elect" in e:
            return True
        elif "qsar" in e:
            return True
        elif "sse" in e:
            return True
        elif "ta" in e:
            return True
        else:
            return False

    def get_data(dfm, metric):
        dfm = dfm.drop("fold", axis=1)
        dfm_count = dfm.apply(np.mean).groupby(by=lambda x: "psekraac" if "lambda" in x or "g-gap" in x else x[:6]).count().to_frame("count")
        dfm_max = dfm.apply(np.mean).groupby(by=lambda x: "psekraac" if "lambda" in x or "g-gap" in x else x[:6]).max().to_frame("max_metric")
        dfm = pd.concat([dfm_max, dfm_count], axis=1)
        dfm["group"] = dfm.index
        dfm["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in dfm["group"]]
        dfm["metric"] = metric
        return dfm

    dfm = pd.concat([get_data(df_f1, "F1"), get_data(df_mcc, "MCC")])

    base = alt.Chart(dfm)

    return alt.layer(
        base.mark_circle(filled=True, size=50, opacity=1.0).encode(
            x=alt.X("group:N", axis=alt.Axis(title=None)),
            y=alt.Y("max_metric", axis=alt.Axis(title=None), scale=alt.Scale(domain=[0.0, 1.0])),
            size=alt.Size("count"),
            color=alt.Color("type:N", scale=alt.Scale(domain=["sequence based", "structure based"],
                                                      range=["#7570b3", "#d95f02"])),
            tooltip="count"
        ),
        base.mark_rule(opacity=0.3).encode(
            x=alt.X("group:N"),
            y="max_metric",
            color="type:N"
        )
    ).facet(
        column=alt.Column("type:N", title=None),
        row=alt.Row("metric:N", title=None)
    ).resolve_scale(
        x='independent'
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