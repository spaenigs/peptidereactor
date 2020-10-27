from more_itertools import chunked
from scipy.spatial import ConvexHull

import pandas as pd
import altair as alt


df_tsne = pd.read_json("data/multiple_datasets/vis/md_tsne/tsne_data.json")

df_convex_hull = df_tsne.groupby(by="dataset") \
    .apply(lambda df: ConvexHull(df.iloc[:, :-1].values)) \
    .to_frame("convex_hull")
df_convex_hull["area"] = df_convex_hull["convex_hull"].apply(lambda h: h.area)
df_convex_hull.sort_values(by="area", inplace=True)

df_convex_hull["ds_field"] = range(1, len(df_convex_hull.index)+1)

range_dict = {}
for ds, s in df_convex_hull.iterrows():
    field_int = s["ds_field"]
    r = [f"{i[0]}-{i[-1]}" for i in chunked(range(1, 51), 12) if field_int in i][0]
    range_dict[ds] = r

df_tsne["range_field"] = df_tsne.dataset.apply(lambda ds: range_dict[ds])

df_tsne["area"] = df_tsne.dataset.apply(lambda ds: df_convex_hull.loc[ds, "area"])

xs, ys = [], []
for i in df_convex_hull.index:
    df_tmp = df_tsne.loc[df_tsne["dataset"].isin([i]), :]
    x = df_tmp.iloc[df_convex_hull.loc[i, "convex_hull"].vertices, 0]
    y = df_tmp.iloc[df_convex_hull.loc[i, "convex_hull"].vertices, 1]
    xs += x.to_list()
    ys += y.to_list()

xs, ys = sorted(xs), sorted(ys)
x_min, x_max, y_min, y_max = xs[0], xs[-1], ys[0], ys[-1]

df_hull = pd.DataFrame()
for ds in df_convex_hull.index:
    hull = df_convex_hull.loc[ds, "convex_hull"]
    df_tmp = df_tsne.loc[df_tsne.dataset == ds, :].copy()
    df_tmp["hull_vertex"] = False
    df_tmp.iloc[hull.vertices, -1] = True
    df_tmp["order"] = -1
    last_idx = len(hull.vertices)+1
    df_tmp.iloc[hull.vertices, -1] = range(1, last_idx)
    tmp = df_tmp.loc[df_tmp.order == 1, :].copy()
    tmp.order += last_idx
    df_hull = pd.concat([df_hull, df_tmp, tmp])

input_dropdown = alt.binding_select(options=df_tsne.range_field.unique())
selection = alt.selection_single(fields=['range_field'], bind=input_dropdown, init={"range_field": "1-12"})

scatterc = alt.Chart().mark_circle(size=3, color="#fdc086").encode(
    x=alt.X(
        "x:Q",
        title="tSNE-1",
        axis=alt.Axis(grid=False, titleFontWeight="normal"),
        scale=alt.Scale(domain=[x_min, x_max])
    ),
    y=alt.Y(
        "y:Q",
        title="tSNE-2",
        axis=alt.Axis(grid=False, titleFontWeight="normal"),
        scale=alt.Scale(domain=[y_min, y_max])
    )
)

print(df_hull.head())

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

alt.layer(scatterc, hullc, data=df_hull).properties(
    width=130,
    height=130
).facet(
    facet=alt.Facet("dataset:N", sort=alt.EncodingSortField("area")),
    columns=4
).transform_filter(
    selection
).add_selection(
    selection
).save("chart.html")




# chart_rows = []
# for c in chunked(range(15), 5):
#
#     print(c)
#
#     chart_column = []
#     for n in df_convex_hull.index[c[0]:c[-1]]:
#         df_tmp = df_tsne.loc[df_tsne["dataset"].isin([n]), :]
#         # url = config["html_dir_out"] + f"tsne_data_{secrets.token_hex(4)}.json"
#         # df_tmp.to_json(url, orient="records")
#         scatterc = alt.Chart(df_tmp).mark_circle(size=3, color="#fdc086").encode(
#             x=alt.X(
#                 "x:Q",
#                 title="tSNE-1",
#                 axis=alt.Axis(grid=False, titleFontWeight="normal"),
#                 scale=alt.Scale(domain=[x_min, x_max])
#             ),
#             y=alt.Y(
#                 "y:Q",
#                 title="tSNE-2",
#                 axis=alt.Axis(grid=False, titleFontWeight="normal"),
#                 scale=alt.Scale(domain=[y_min, y_max])
#             )
#         ).properties(
#             height=130,
#             width=130
#         )
#
#         hull = df_convex_hull.loc[n, "convex_hull"]
#         points = df_tmp.iloc[hull.vertices, :]
#         source = pd.DataFrame({"x": points["x"], "y": points["y"]})
#         source = pd.concat([source, source.iloc[:1, :]])
#         source["order"] = range(1, source.shape[0] + 1)
#         source.index = range(0, source.shape[0])
#         # url_hull = config["html_dir_out"] + f"tsne_fhull_{secrets.token_hex(4)}.json"
#         # source.to_json(url_hull, orient="records")
#         hullc = alt.Chart(source).mark_line(
#             color="#386cb0",
#             strokeDash=[5, 3],
#             strokeWidth=1
#         ).encode(
#             x="x:Q",
#             y="y:Q",
#             order="order:O"
#         )
#
#         df_text = pd.DataFrame({"x": [0], "y": [0], "area": [hull.area]})
#         textc = alt.Chart(df_text).mark_text().encode(
#             x="x:Q",
#             y="y:Q",
#             text="text:N"
#         ).transform_calculate(
#             text="join(['area=', round(datum.area)], '')"
#         )
#
#         chart_column += [alt.layer(scatterc, hullc, textc, title=alt.TitleParams(text=n, anchor="middle"))]
#
#     if len(chart_column) != 0:
#         chart_rows += [alt.hconcat(*chart_column)]
#
# alt.vconcat(*chart_rows).save("chart.html")
