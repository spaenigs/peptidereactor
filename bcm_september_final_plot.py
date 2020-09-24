import pandas as pd
import altair as alt

source = pd.read_csv("source.csv")

sorted_names = source.groupby("seq_name").apply(lambda df: df.iloc[0, 5]).sort_values(ascending=False).index
d = dict((k,v) for k,v in zip(sorted_names, range(len(sorted_names))))

source["rank"] = source["seq_name"].replace(d)
source.sort_values(["rank", "y"], inplace=True)

xvals = source["x"].unique()
seq_names = source["seq_name"].unique()

input_dropdown = alt.binding_select(options=list(seq_names) + [None])
selection = alt.selection_single(fields=['seq_name'], bind=input_dropdown, name="_")

alt.Chart(source).mark_point(filled=True, shape="square", size=300, opacity=1.0).encode(
    x=alt.X("x:N", title=None, sort=alt.Sort(alt.SortArray(xvals))),
    y=alt.Y("y:N", title="Lab results vs. prediction", sort=alt.Sort(alt.SortArray(source.sort_values(["rank", "y"])["y"].unique()))),
    color=alt.Color("coupl:N", title="% Ac"),
    tooltip="seq:N"
).add_selection(
    selection
).transform_filter(
    selection
).configure_axis(
    labelLimit=350,
    labelFontSize=12,
    titlePadding=150
).properties(
    title=["Comparison of difficult couplings observed in the lab vs. predictions",
           "Window length: 7, cut-off: 1.5 %"],
    width=700
).save("chart.html")