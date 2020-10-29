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

ds = alt.Chart().mark_text(fontSize=12, lineBreak="\n").encode(
    text="dataset:N"
).properties(
    height=80,
    width=125
)

ss = alt.Chart().mark_text(fontSize=12).encode(
    text="seq_size:O",
).properties(
    height=80,
    width=125
)

t = alt.Chart().mark_text(
    dy=-30,
    fontSize=12,
    lineBreak="\n",
).encode(
    text="desc:N"
)

ref = alt.Chart().transform_calculate(
    url="" + alt.datum.ref
).mark_text(
    dy=-30,
    fontSize=12,
    lineBreak="\n"
).encode(
    text="ref:N",
    href="url:N",
    tooltip="url:N"
)

table = alt.hconcat(ds, ss, t, ref, data=snakemake.input[0]).transform_filter(
    selection
)

scatter = alt.Chart(
    snakemake.input[0],
    title="Multiple datasets"
).mark_point(filled=True, size=60).encode(
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
).add_selection(
    selection
).properties(
    height=250,
    width=250
).interactive()

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
).transform_filter(
    selection
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
).transform_filter(
    selection
)

textc = alt.Chart().mark_text().encode(
    x="x:Q",
    y="y:Q",
    text="text:N"
).transform_calculate(
    text="join(['area=', round(datum.area)], '')"
).transform_filter(
    (alt.datum.x == 0) and (alt.datum.y == 0)
).transform_filter(
    selection
)

tsnec = alt.layer(
    scatterc, hullc, textc,
    data=snakemake.input[1],
    title="Single dataset"
).properties(
    height=250,
    width=250
)

chart_json = alt.vconcat(alt.hconcat(scatter, tsnec), table).to_json(indent=None)

with open(snakemake.output[0], "w") as f:
    f.write(chart_json)
    f.flush()