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
df_sb.Encoding = ""
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
names_seq = names_seq[1:] + [names_seq[0]]
names_str = list(source.loc[source.type == "structure based", "Encoding"].sort_values().unique())

x_axis_config = alt.Axis(labelAngle=-45)
y_axis_config = alt.Axis(tickCount=len(source.Encoding))

sort_x_axis = alt.Sort(alt.SortArray(names_seq + names_str))
sort_y_axis = alt.Sort(alt.SortArray(names_imbalanced))

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

chart1 = alt.Chart(url).mark_rect(
    size=RECT_SIZE,
    stroke='black',
    strokeWidth=0.2
).encode(
    x=alt.X("Encoding:N", axis=x_axis_config, sort=sort_x_axis),
    y=alt.Y("Dataset:N", axis=y_axis_config, sort=sort_y_axis),
    color=color,
    tooltip=tooltip
)

chart2 = alt.Chart(url).mark_rect(size=RECT_SIZE).encode(
    x=alt.X('Encoding:N', axis=x_axis_config, sort=sort_x_axis),
    y=alt.Y('Dataset:N', axis=y_axis_config, sort=sort_y_axis),
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

sep_chart = alt.Chart(url).mark_rect(strokeOpacity=1.0).encode(
    x=alt.X('Encoding:N', axis=x_axis_config, sort=sort_x_axis),
    y=alt.Y('Dataset:N', sort=sort_y_axis),
    color=alt.value("#f0f0f0")
)

chart3 = sep_chart.transform_filter(
    alt.datum.Dataset == ""
)

chart4 = sep_chart.transform_filter(
    alt.datum.Encoding == ""
)

chart = alt.layer(
    chart1, chart2, chart3, chart4
).resolve_scale(
    color="independent"
)

joblib.dump(chart, snakemake.output[0])

