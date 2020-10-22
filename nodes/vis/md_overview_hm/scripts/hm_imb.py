import pandas as pd
import altair as alt

import joblib

RECT_SIZE = 13

source = pd.read_json(snakemake.input[0])

df_ds = source.loc[source.Dataset == "cpp_sanders", :].copy(deep=True)
df_ds.Dataset = ""
df_ds.F1 = "separator"
df_ds.is_imbalanced = 0.35

df_sb = source.loc[source.Encoding == "zscale", :].copy(deep=True)
df_sb.Encoding = "zzz"
df_sb.F1 = "separator"

source = pd.concat([source, df_ds, df_sb]).sort_values(by="is_imbalanced")

url = snakemake.input[0].replace("source", "imb")
source.to_json(url, orient="records")

names_imbalanced = list(reversed(source["Dataset"].unique()))
names_seq = list(source.loc[source.type == "sequence based", "Encoding"].sort_values().unique())
names_str = list(source.loc[source.type == "structure based", "Encoding"].sort_values().unique())

axis = alt.Axis(
    tickCount=len(source.Encoding),
    labelExpr="datum.label == 'zzz' ? null : datum.label"
)

sort_x_axis = alt.Sort(alt.SortArray(names_seq + names_str))
sort_y_axis = alt.Sort(alt.SortArray(names_imbalanced))

tooltip = ["Encoding:N", "Dataset:N", "F1:Q", "is_imbalanced:Q"]

chart1 = alt.Chart(url).mark_rect().encode(
    x=alt.X('Encoding:N', axis=axis, sort=sort_x_axis),
    y=alt.Y('Dataset:N', sort=sort_y_axis),
    color=alt.Color('F1:Q', scale=alt.Scale(
        domain=[0.0, 1.0], range=["#a6bddb", "#023858"]
    )),
    tooltip=tooltip
)

chart2 = alt.Chart(url).mark_rect(size=RECT_SIZE).encode(
    x=alt.X('Encoding:N', axis=axis, sort=sort_x_axis),
    y=alt.Y('Dataset:N', sort=sort_y_axis),
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

chart = (chart1 + chart2)

joblib.dump(chart, snakemake.output[0])
