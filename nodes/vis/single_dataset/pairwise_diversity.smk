import pandas as pd
import altair as alt
import numpy as np

from nodes.vis.single_dataset.scripts.utils import *


def get_crap_combination():
    i_crap, c_crap, div_crap = \
             "cgr_res_20_sf_0.5", "electrostatic_hull_6", df_div.loc["cgr_res_20_sf_0.5", "electrostatic_hull_6"]
         df1_crap, df2_crap, df1_crap_class = get_data(i_crap, c_crap)
         df_crap = ravel_and_annotate(df1_crap, df2_crap, df1_crap_class,
                                     cat=f"(a) crap", div=div_crap, e1=i_crap, e2=c_crap,
                                     f1_e1=medians[i_crap], f1_e2=medians[c_crap])

TOKEN = config["token"]

rule pairwise_div_transform_data:
    input:
         config["metrics_dir_in"] + "f1.csv",
         config["similarity_dir_in"] + "{comparision}/diversity.csv"
    output:
         f"data/temp/{TOKEN}/{{comparision}}.csv"
    run:
         df_f1 = pd.read_csv(input[0], index_col=0)
         df_div = pd.read_csv(input[1], index_col=0)

         medians = df_f1.apply(np.median).sort_values(ascending=False)
         indices = medians.index.tolist()

         if wildcards.comparision == "seq_vs_str":
             indices_1 = indices[:15]
             indices_2 = [i for i in indices if is_struc_based(i)][:5]
         else:
             indices_1 = indices[:15]
             indices_2 = indices[:15]

         df_div_sub = df_div.loc[indices_1, indices_2]

         sorted([(s, i, df_div.loc[s, i]) for i, s in df_div.idxmin().items()], key=lambda x: x[2])





