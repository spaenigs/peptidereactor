from more_itertools import chunked
from iFeature import AAC
from scipy.spatial import ConvexHull
from modlamp.core import read_fasta
from sklearn.manifold import TSNE

import pandas as pd
import altair as alt

import re
import secrets

TOKEN = config["token"]

rule transform_tsne_data:
    input:
         fastas=config["fastas_in"],
         classes=config["classes_in"]
    output:
         config["html_dir_out"] + f"data/tsne/data.json"
    run:
        df_res = pd.DataFrame()
        for fasta_path, class_path in zip(input.fastas, input.classes):
            seqs, names = read_fasta(fasta_path)
            with open(class_path) as f:
                classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
            seq_tuples = [[name, tup[0]]
                          for name, tup in zip(names, zip(seqs, classes))
                          if tup[1] == 1]
            df_tmp = pd.DataFrame([res[1:] for res in AAC.AAC(seq_tuples, order=None)][1:])
            df_tmp["dataset"] = re.findall("data/(.*?)/", fasta_path)[0]
            df_res = pd.concat([df_res, df_tmp])

        # df_res.columns = [res[1:] for res in AAC.AAC(seq_tuples, order=None)][0] + [df_res.columns[-1]]
        X_embedded = TSNE(n_components=2, n_jobs=10).fit_transform(df_res.iloc[:, :-1].values)

        df_tsne = pd.DataFrame(X_embedded)
        df_tsne.columns = ["x", "y"]
        df_tsne["dataset"] = df_res["dataset"].to_list()

        df_tsne.to_json(output[0], orient="records")

rule tsne_chart:
    input:
         config["html_dir_out"] + f"data/tsne/data.json"
    output:
         temp(f"data/temp/{TOKEN}/tsne.json")
    run:
         df_tsne = pd.read_json(input[0])
         url = input[0].replace(config["html_dir_out"], "")

         df_convex_hull = df_tsne.groupby(by="dataset")\
            .apply(lambda df: ConvexHull(df.iloc[:, :-1].values))\
            .to_frame("convex_hull")
         df_convex_hull["area"] = df_convex_hull["convex_hull"].apply(lambda h: h.area)
         df_convex_hull.sort_values(by="area", inplace=True)

         xs, ys = [], []
         for i in df_convex_hull.index:
             df_tmp = df_tsne.loc[df_tsne["dataset"].isin([i]), :]
             x = df_tmp.iloc[df_convex_hull.loc[i, "convex_hull"].vertices, 0]
             y = df_tmp.iloc[df_convex_hull.loc[i, "convex_hull"].vertices, 1]
             xs += x.to_list()
             ys += y.to_list()

         xs, ys = sorted(xs), sorted(ys)
         x_min, x_max,y_min, y_max = xs[0], xs[-1], ys[0], ys[-1]

         chart_rows = []
         for c in chunked(range(15), 5):

             chart_column = []
             for n in df_convex_hull.index[c[0]:c[-1]]:

                 df_tmp = df_tsne.loc[df_tsne["dataset"].isin([n]), :]
                 filename = config["html_dir_out"] + f"data/tsne/data_{secrets.token_hex(4)}.json"
                 df_tmp.to_json(filename, orient="records")
                 url = filename.replace(config["html_dir_out"], "")
                 scatterc = alt.Chart(url).mark_circle(size=3, color="#fdc086").encode(
                     x=alt.X(
                         "x:Q",
                         title="tSNE-1",
                         axis=alt.Axis(grid=False),
                         scale=alt.Scale(domain=[x_min, x_max])
                     ),
                     y=alt.Y(
                         "y:Q",
                         title="tSNE-2",
                         axis=alt.Axis(grid=False),
                         scale=alt.Scale(domain=[y_min, y_max])
                     ),
                     tooltip="dataset:N"
                  ).properties(
                     height=130,
                     width=130
                 )

                 hull = df_convex_hull.loc[n, "convex_hull"]
                 points = df_tmp.iloc[hull.vertices, :]
                 source = pd.DataFrame({"x": points["x"], "y": points["y"]})
                 source = pd.concat([source, source.iloc[:1, :]])
                 source["order"] = range(1, source.shape[0]+1)
                 source.index = range(0, source.shape[0])
                 filename = config["html_dir_out"] + f"data/tsne/hull_{secrets.token_hex(4)}.json"
                 source.to_json(filename, orient="records")
                 url_hull = filename.replace(config["html_dir_out"], "")
                 hullc = alt.Chart(url_hull).mark_line(
                     color="#386cb0",
                     strokeDash=[5, 3],
                     strokeWidth=1
                 ).encode(
                     x="x:Q",
                     y="y:Q",
                     order="order:O"
                 )

                 df_text = pd.DataFrame({"x": [0], "y": [0], "area": [hull.area]})
                 textc = alt.Chart(df_text).mark_text().encode(
                     x="x:Q",
                     y="y:Q",
                     text="text:N"
                 ).transform_calculate(
                     text="join(['area=', round(datum.area)], '')"
                 )

                 chart_column += [alt.layer(scatterc, hullc, textc, title=alt.TitleParams(text=n, anchor="middle"))]

             chart_rows += [alt.hconcat(*chart_column)]

         chart_json = alt.vconcat(*chart_rows).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()

