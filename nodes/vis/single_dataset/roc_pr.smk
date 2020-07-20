from functools import reduce

import altair as alt

import re
import joblib

from nodes.vis.single_dataset.scripts.roc_pr_curve import *
from nodes.vis.single_dataset.scripts.utils import *


def transform_curve_data(encodings, dataset, curve):
    test_dfs, prob_dfs, names = [], [], []
    for e in encodings:
        test_dfs += [pd.read_csv(f"data/{dataset}/benchmark/single/y_true_cv_{e}.csv", index_col=0)]
        prob_dfs += [pd.read_csv(f"data/{dataset}/benchmark/single/y_prob_cv_{e}.csv", index_col=0)]
        names += [e]
    df = reduce(
        lambda a, b: pd.concat([a, b]),
        map(get_roc_data if curve == "roc" else get_pr_data, test_dfs, prob_dfs, names)
    )
    df["type"] = ["structure based" if is_struc_based(e) else "sequence based"
                  for e in df["Encoding"]]
    return df


TOKEN = config["token"]

rule metric_curves_transform_data:
    input:
         config["metrics_dir_in"] + "f1.csv"
    output:
         temp(f"data/temp/{TOKEN}/{{curve}}_data_{{topencs}}.csv")
    run:
         df_f1 = pd.read_csv(input[0], index_col=0)

         indices = df_f1.apply(np.median)\
             .sort_values(ascending=False)\
             .index.tolist()

         dataset = re.findall("data/(.*?)/", config["html_out"])[0]

         if wildcards.topencs == "t6":
            top_six_encodings = indices[:6]
            df = transform_curve_data(top_six_encodings, dataset, wildcards.curve)
            df.to_csv(output[0])
         else:
            top_three_seq = indices[:3]
            top_three_str = [i for i in indices if is_struc_based(i)][:3]
            df = transform_curve_data(top_three_seq + top_three_str, dataset, wildcards.curve)

         df.to_csv(output[0])

rule create_sub_chart:
    input:
         f"data/temp/{TOKEN}/{{curve}}_data_{{topencs}}.csv"
    output:
         temp(f"data/temp/{TOKEN}/{{curve}}_{{topencs}}.joblib")
    run:
         random_guess_line = alt.Chart(
                 pd.DataFrame({"x": [0.0, 1.0],
                               "y": [0.0, 1.0]})
             ).mark_line(
                 color="lightgrey",
                 strokeDash=[3, 1]
             ).encode(
                 x="x",
                 y="y"
             )

         path = input[0]
         df = pd.read_csv(path, index_col=0)

         if "roc" in path and "t6" in path:
             chart = alt.Chart(df).mark_line().encode(
                 x=alt.X("x:Q", axis=alt.Axis(title=None)),
                 y=alt.Y("y:Q", axis=alt.Axis(title="Sensitivity")),
                 color=alt.Color("Encoding:N")
             ) + random_guess_line
         elif "pr" in path and "t6" in path:
             chart = alt.Chart(df).mark_line().encode(
                 x=alt.X("x:Q", axis=alt.Axis(title=None)),
                 y=alt.Y("y:Q", axis=alt.Axis(title="Precision")),
                 color="Encoding:N"
             )
         elif "roc" in path and "t33" in path:
             chart = alt.Chart(df).mark_line().encode(
                 x=alt.X("x:Q", axis=alt.Axis(title="1 - Specificity")),
                 y=alt.Y("y:Q", axis=alt.Axis(title="Sensitivity")),
                 color="Encoding:N",
                 strokeDash="type:N"
             ) + random_guess_line
         else:
             chart = alt.Chart(df).mark_line().encode(
                 x=alt.X("x:Q", axis=alt.Axis(title="Recall")),
                 y=alt.Y("y:Q", axis=alt.Axis(title="Precision")),
                 color="Encoding:N",
                 strokeDash="type:N"
             )

         joblib.dump(chart, output[0])

rule create_roc_pr_chart:
    input:
         expand(f"data/temp/{TOKEN}/{{curve}}_{{topencs}}.joblib",
                curve=["roc", "pr"], topencs=["t6", "t33"])
    output:
         temp(f"data/temp/{TOKEN}/roc_pr.json")
    run:
         paths = list(input)
         c1 = joblib.load(path(paths, "roc", "t6"))
         c2 = joblib.load(path(paths, "pr", "t6"))
         cA = joblib.load(path(paths, "roc", "t33"))
         cB = joblib.load(path(paths, "pr", "t33"))

         chart = (c1 | c2) & (cA | cB)

         with open(output[0], "w") as f:
             f.write(chart.to_json(indent=None))
             f.flush()

