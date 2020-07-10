import pandas as pd
import altair as alt
import numpy as np


def overview_chart(df_f1, df_mcc):

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
        dfm_count = dfm.apply(np.mean).groupby(
            by=lambda x: "psekraac" if "lambda" in x or "g-gap" in x else x[:6]).count().to_frame("count")
        dfm_max = dfm.apply(np.mean).groupby(
            by=lambda x: "psekraac" if "lambda" in x or "g-gap" in x else x[:6]).max().to_frame("max_metric")
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