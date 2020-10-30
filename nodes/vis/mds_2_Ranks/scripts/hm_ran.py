import pandas as pd
import altair as alt

import joblib


def annotate(df):
    df_tmp = df.sort_values("F1", ascending=False)
    df_tmp["rank"] = [str(i) if i < 4 else ">3" for i, e in enumerate(df_tmp["F1"], start=1)]
    return df_tmp


RECT_SIZE = 13

source = pd.read_json(snakemake.input[0])
source = source.groupby("Dataset").apply(annotate).reset_index(drop=True)

df_sb = source.loc[source.Encoding == "zscale", :].copy(deep=True)
df_sb.Encoding = "zzz"
df_sb.F1 = "separator"
df_sb["rank"] = ">3"

source = pd.concat([source, df_sb])

df_ds = source.loc[source.Dataset == "cpp_sanders", :].copy(deep=True)
df_ds.Dataset = ""
df_ds.F1 = "separator"
df_ds.is_imbalanced = 0.35
df_ds["rank"] = ">3"

source = pd.concat([source, df_ds]).sort_values(by="is_imbalanced")

names_imbalanced = list(reversed(source["Dataset"].unique()))
names_seq = list(source.loc[source.type == "sequence based", "Encoding"].sort_values().unique())
names_str = list(source.loc[source.type == "structure based", "Encoding"].sort_values().unique())

axis = alt.Axis(
    tickCount=len(source.Encoding),
    labelExpr="datum.label == 'zzz' ? null : datum.label"
)

sort_y_axis = alt.Sort(alt.SortArray(names_seq + names_str))
sort_x_axis = alt.Sort(alt.SortArray(names_imbalanced))

color = alt.condition(
    alt.datum.rank == 4,
    alt.value("white"),
    alt.Color(
        "rank:N",
        title="Rank",
        scale=alt.Scale(
            domain=["1", "2", "3", ">3"],
            range=["#023858", "#0570b0", "#74a9cf", "white"]
        ),
        legend=alt.Legend(symbolStrokeWidth=0.2)
    )
)

tooltip = ["Encoding:N", "Dataset:N", "F1:Q", "is_imbalanced:Q"]

url = snakemake.input[0].replace("source", "ran")
source.to_json(url, orient="records")

chart1 = alt.Chart(url).mark_rect(stroke='black', strokeWidth=0.2).encode(
    y=alt.Y('Encoding:N', axis=axis, sort=sort_y_axis),
    x=alt.X('Dataset:N', axis=alt.Axis(labelAngle=-45), sort=sort_x_axis),
    color=color,
    tooltip=tooltip
)

chart2 = alt.Chart(url).mark_rect(size=RECT_SIZE).encode(
    y=alt.Y('Encoding:N', axis=axis, sort=sort_y_axis),
    x=alt.X('Dataset:N', axis=alt.Axis(labelAngle=-45), sort=sort_x_axis),
    color=alt.Color(
        'F1_new:N',
        title="Value",
        scale=alt.Scale(domain=["NA"], range=["#a6611a"])
    )
).transform_calculate(
    F1_new="datum.F1 == null ? 'NA' : 'NA'"
).transform_filter(
    alt.datum.F1 == None
).properties(
    height={"step": RECT_SIZE},
    width={"step": RECT_SIZE}
)

sep_chart = alt.Chart(url).mark_rect().encode(
    y=alt.Y('Encoding:N', axis=axis, sort=sort_y_axis),
    x=alt.X('Dataset:N', sort=sort_x_axis),
    color=alt.value("#f0f0f0")
)

chart3 = sep_chart.transform_filter(
    alt.datum.Dataset == ""
)

chart4 = sep_chart.transform_filter(
    alt.datum.Encoding == "zzz"
)

chart = (chart1 + chart2 + chart3 + chart4).resolve_scale(color="independent")

joblib.dump(chart, snakemake.output[0])

