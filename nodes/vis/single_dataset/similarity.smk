import altair as alt
import pandas as pd

import joblib
import os

from nodes.vis.single_dataset.scripts.utils import *

TOKEN = config["token"]

SIMILARITY_DIR_IN = os.path.commonpath(
    [config["similarity_dir_group_1_in"], config["similarity_dir_group_2_in"]]) + "/"

rule similarity_transform_data:
    input:
         SIMILARITY_DIR_IN + "{comparision}/{metric}.csv"
    output:
         temp(f"data/temp/{TOKEN}/{{comparision}}_{{metric}}.data")
    run:
         df = pd.read_csv(input[0], index_col=0)

         row_indices = cluster(df.values, axis=0)
         col_indices = cluster(df.values, axis=1)
         heatmap_data = df.iloc[row_indices, col_indices]

         x, y = np.meshgrid(range(0, heatmap_data.shape[1]), range(0, heatmap_data.shape[0]))
         source = pd.DataFrame({"x": x.ravel(), "y": y.ravel(), "Similarity": heatmap_data.values.ravel()})

         source["Encoding1"] = source["x"].apply(lambda i: heatmap_data.columns[i])
         source["Encoding2"] = source["y"].apply(lambda i: heatmap_data.index[i])

         source["Similarity_cat"] = source.Similarity.apply(
             lambda x: "<= 0.1" if x <= 0.1 else "<= 0.3" if x <= 0.3 else "<= 0.6" if x <= 0.6 else "<= 1.0")

         source.to_csv(output[0])

rule create_heatmap:
    input:
         f"data/temp/{TOKEN}/{{comparision}}_{{metric}}.data"
    output:
         temp(f"data/temp/{TOKEN}/{{comparision}}_{{metric}}.hmjl")
    run:
         d, r = ["<= 0.1", "<= 0.3", "<= 0.6", "<= 1.0"], ["white", "gainsboro", "grey", "black"]

         x_config = alt.Axis(labels=False, ticks=False)
         y_config = alt.Axis(labels=False, ticks=False)

         comp, met = wildcards.comparision, wildcards.metric

         show_x_title, show_y_title = True, True
         if comp == "all_vs_all" and met == "diversity":
             show_x_title = False
         elif comp == "all_vs_all" and met == "phi":
             show_x_title, show_y_title = True, True
         elif comp == "seq_vs_str" and met == "phi":
             show_y_title = True
         else:
             show_x_title, show_y_title = True, True

         if not show_x_title:
             x_config = alt.Axis(labels=False, ticks=False, title=None)
         if not show_y_title:
             y_config = alt.Axis(labels=False, ticks=False, title=None)

         source = pd.read_csv(input[0], index_col=0)

         chart = alt.Chart(source).mark_rect().encode(
             x=alt.X('x:O', title="Encoding 2", axis=x_config),
             y=alt.Y('y:O', title="Encoding 1", axis=y_config),
             color=alt.Color(
                 "Similarity_cat:N",
                 scale=alt.Scale(domain=d, range=r),
                 legend=alt.Legend(title="Similarity Range")
             ),
             tooltip=["Encoding1", "Encoding2", "Similarity"]
         ).properties(
             width=600,
             height=600
         )

         joblib.dump(chart, output[0])

rule create_similarity_chart:
    input:
         expand(f"data/temp/{TOKEN}/{{comparision}}_{{metric}}.hmjl",
                comparision=["all_vs_all", "seq_vs_str"],
                metric=["diversity", "phi"])
    output:
         temp(f"data/temp/{TOKEN}/similarity.json")
    run:
         paths = list(input)

         hm1 = joblib.load(path(paths, "all_vs_all", "diversity"))
         hm2 = joblib.load(path(paths, "all_vs_all", "phi"))
         hm3 = joblib.load(path(paths, "seq_vs_str", "diversity"))
         hm4 = joblib.load(path(paths, "seq_vs_str", "phi"))

         t1 = alt.hconcat(hm1, hm2, title=alt.TitleParams(text="All vs. all", anchor="middle"))
         t2 = alt.hconcat(hm3, hm4, title=alt.TitleParams(text="Seq. vs. struc.", anchor="middle"))

         alt.data_transformers.disable_max_rows()
         chart_json = (t1 & t2).configure_view(strokeWidth=0).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()