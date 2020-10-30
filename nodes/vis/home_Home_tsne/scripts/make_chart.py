import altair as alt

import yaml


with open(snakemake.input[2]) as f:
    x_min, x_max, y_min, y_max = yaml.safe_load(f)

selection = alt.selection_single(
    fields=["dataset"],
    init={"dataset": "hiv_ddi"},
    empty="none"
)

axis = alt.Axis(grid=False, titleFontWeight="normal")

### chart1

scatter = alt.Chart(
    snakemake.input[0], title="Multiple datasets"
).mark_point(
    filled=True, size=60
).encode(
    x=alt.X(
        "hours:Q",
        title="log(Computation time) (h)",
        scale=alt.Scale(type='log'),
        axis=axis
    ),
    y=alt.Y(
        "seq_size:Q",
        title="log(# of sequences)",
        sort="-x",
        scale=alt.Scale(type='log'),
        axis=axis
    ),
    tooltip="dataset:N",
    color=alt.condition(selection, alt.value("#4C78A8"), alt.value("lightgrey"))
).interactive()

### chart2

scatterc = alt.Chart().mark_circle(size=10, color="#fdc086").encode(
    x=alt.X(
        "x:Q",
        title="tSNE-1",
        axis=axis,
        scale=alt.Scale(domain=[x_min, x_max])
    ),
    y=alt.Y(
        "y:Q",
        title="tSNE-2",
        axis=axis,
        scale=alt.Scale(domain=[y_min, y_max])
    )
)

hullc = alt.Chart().mark_line(
    color="#386cb0",
    strokeDash=[5, 3],
    strokeWidth=1
).encode(
    x="x:Q",
    y="y:Q",
    order="order:O"
).transform_filter(
    alt.datum.hull_vertex == True
)

textc = alt.Chart().mark_text().encode(
    x="x:Q",
    y="y:Q",
    text="text:N"
).transform_calculate(
    text="join(['area=', round(datum.area)], '')"
).transform_filter(
    (alt.datum.x == 0) and (alt.datum.y == 0)
)

tsnec = alt.layer(
    scatterc, hullc, textc,
    data=snakemake.input[1],
    title="Single dataset"
)

### chart3

ds = alt.Chart().mark_text(fontSize=12, lineBreak="\n").encode(
    text="dataset:N"
)

ss = alt.Chart().mark_text(fontSize=12).encode(
    text="seq_size:O",
)

t = alt.Chart().mark_text(
    dy=-50,
    fontSize=12,
    lineBreak="\n",
).encode(
    text="desc:N"
)

ref = alt.Chart().transform_calculate(
    url="" + alt.datum.ref
).mark_text(
    dy=-50,
    fontSize=12,
    lineBreak="\n"
).encode(
    text="ref:N",
    href="url:N",
    tooltip="url:N"
).properties(
    height=120,
    width=125
)

# chart

chart1 = scatter.add_selection(
    selection
).properties(
    height=250,
    width=250
)

chart2 = alt.layer(
    scatterc, hullc, textc,
    data=snakemake.input[1],
    title="Single dataset"
).transform_filter(
    selection
).properties(
    height=250,
    width=250
)

chart3 = alt.hconcat(
    ds, ss, t, ref, data=snakemake.input[0]
).transform_filter(
   selection
)

chart = alt.vconcat(
    alt.hconcat(chart1, chart2), chart3
).resolve_scale(
    y="independent"
)

with open(snakemake.output[0], "w") as f:
    f.write(chart.to_json(indent=None))
    f.flush()