import altair as alt

import joblib
import os

from nodes.vis.single_dataset.scripts.utils import *
from nodes.vis.single_dataset.scripts.pairwise_diversity import *

TOKEN = config["token"]

SIMILARITY_DIR_IN = os.path.commonpath(
    [config["similarity_dir_group_1_in"], config["similarity_dir_group_2_in"]]) + "/"

ENSEMBLE_CV_DIR_IN = os.path.commonpath(
    [config["ensemble_cv_group_1a_in"], config["ensemble_cv_group_2a_in"],
     config["ensemble_cv_group_1b_in"], config["ensemble_cv_group_2b_in"]]) + "/"

rule pairwise_div_transform_data:
    input:
         config["metrics_dir_in"] + "f1.csv",
         SIMILARITY_DIR_IN + "{comparision}/diversity.csv",
         ENSEMBLE_CV_DIR_IN + "{comparision}/"
    output:
         temp(f"data/temp/{TOKEN}/{{comparision}}.csv")
    run:
         df_f1 = pd.read_csv(input[0], index_col=0)
         df_div = pd.read_csv(input[1], index_col=0)

         medians = df_f1.apply(np.median).sort_values(ascending=False)
         indices = medians.index.tolist()

         if wildcards.comparision == "seq_vs_str":
             indices_1 = [i for i in indices if not is_struc_based(i)][:15]
             indices_2 = [i for i in indices if is_struc_based(i)][:5]
         else:
             indices_1 = indices[:15]
             indices_2 = indices[15:30]

         df_div_sub = df_div.loc[indices_1, indices_2]

         df_crap = get_crap_combination(input[2], df_div, medians)
         df_low = get_low_combination(input[2], df_div_sub, medians)
         df_mid = get_mid_combination(input[2], df_div_sub, medians)
         df_high = get_high_combination(input[2], df_div_sub, medians)
         df_highest_f1 = \
             get_highest_f1_combination(input[2], indices_1[0], indices_2[0], df_div_sub, medians)

         df = pd.concat([df_crap, df_low, df_mid, df_high, df_highest_f1])
         df.dropna().to_csv(output[0])

rule create_scatter_chart:
    input:
         f"data/temp/{TOKEN}/{{comparision}}.csv"
    output:
         temp(f"data/temp/{TOKEN}/{{comparision}}.scjl")
    run:
         df = pd.read_csv(input[0], index_col=0)

         if wildcards.comparision == "all_vs_all":
             title = "All vs. all encodings"
         else:
             title = "Sequence vs. structure based encoding"

         chart = alt.Chart(df).mark_point(text=title, filled=True, size=5, opacity=1).encode(
             x=alt.X("x", title=f"Predicted probability e1", scale=alt.Scale(domain=[0.0, 1.0])),
             y=alt.Y("y", title=f"Predicted probability e2", scale=alt.Scale(domain=[0.0, 1.0])),
             color=alt.Color(
                 "class:N",
                 scale=alt.Scale(domain=[0, 1], range=["#7570b3", "#d95f02"])),
             column=alt.Column(
                 "diversity",
                 header=alt.Header(labelOrient='left'),
                 title=None
             )
         ).properties(
             width=250,
             height=250
         )

         chart = alt.hconcat(chart, title=alt.TitleParams(
             text=title + " (increasing diversity w.r.t. performance)",
             anchor="middle"
         ))

         joblib.dump(chart, output[0])

rule create_pairwise_diversity_chart:
    input:
         expand(f"data/temp/{TOKEN}/{{comparision}}.scjl",
                comparision=["all_vs_all", "seq_vs_str"])
    output:
         temp(f"data/temp/{TOKEN}/pairwise_diversity.json")
    run:
         paths = list(input)

         sc1 = joblib.load(path(paths, "seq_vs_str", ""))
         sc2 = joblib.load(path(paths, "all_vs_all", ""))

         alt.data_transformers.disable_max_rows()

         chart_json = (sc1 & sc2).to_json(indent=None)

         with open(output[0], "w") as f:
             f.write(chart_json)
             f.flush()
