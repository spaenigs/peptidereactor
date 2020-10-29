from more_itertools import chunked
from scipy.spatial import ConvexHull

import pandas as pd

import yaml


def create_hull_data(df_tsne):
    df_convex_hull = df_tsne.groupby(by="dataset") \
        .apply(lambda df: ConvexHull(df.iloc[:, :-1].values)) \
        .to_frame("convex_hull")
    df_convex_hull["area"] = df_convex_hull["convex_hull"].apply(lambda h: h.area)
    df_convex_hull.sort_values(by="area", inplace=True)
    df_convex_hull["ds_field"] = range(1, len(df_convex_hull.index) + 1)
    return df_convex_hull


def get_range_dict(df_convex_hull):
    range_dict = {}
    for ds, s in df_convex_hull.iterrows():
        field_int = s["ds_field"]
        r = [f"{i[0]}-{i[-1]}" for i in chunked(range(1, 51), 12) if field_int in i][0]
        range_dict[ds] = r
    return range_dict


def get_axis_range(df_convex_hull):
    xs, ys = [], []
    for i in df_convex_hull.index:
        df_tmp = df_tsne.loc[df_tsne["dataset"].isin([i]), :]
        x = df_tmp.iloc[df_convex_hull.loc[i, "convex_hull"].vertices, 0]
        y = df_tmp.iloc[df_convex_hull.loc[i, "convex_hull"].vertices, 1]
        xs += x.to_list()
        ys += y.to_list()
    xs, ys = sorted(xs), sorted(ys)
    x_min, x_max, y_min, y_max = xs[0], xs[-1], ys[0], ys[-1]
    return x_min, x_max, y_min, y_max


def combine_data(df_tsne, df_convex_hull):
    df_hull = pd.DataFrame()
    for ds in df_convex_hull.index:
        hull = df_convex_hull.loc[ds, "convex_hull"]
        df_tmp = df_tsne.loc[df_tsne.dataset == ds, :].copy()
        df_tmp["hull_vertex"] = False
        df_tmp.iloc[hull.vertices, -1] = True
        df_tmp["order"] = -1
        last_idx = len(hull.vertices) + 1
        df_tmp.iloc[hull.vertices, -1] = range(1, last_idx)
        tmp = df_tmp.loc[df_tmp.order == 1, :].copy()
        tmp.order += last_idx
        df_text = tmp.copy()
        df_text.x, df_text.y, df_text.hull_vertex = 0, 0, False
        df_hull = pd.concat([df_hull, df_tmp, tmp, df_text])
    df_hull.reset_index(inplace=True)
    return df_hull


df_tsne = pd.read_csv(snakemake.input[0], index_col=0)

df_convex_hull = create_hull_data(df_tsne)

range_dict = get_range_dict(df_convex_hull)
df_tsne["range_batch"] = df_tsne.dataset.apply(lambda ds: range_dict[ds])
df_tsne["area"] = df_tsne.dataset.apply(lambda ds: df_convex_hull.loc[ds, "area"])

df_hull = combine_data(df_tsne, df_convex_hull)
df_hull.to_json(snakemake.output[0], orient="records")

with open(snakemake.output[1], "w") as f:
    yaml.safe_dump(get_axis_range(df_convex_hull), f)
