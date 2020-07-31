import pandas as pd
import numpy as np
import altair as alt

from glob import glob

# stack = np.dstack((df_seq.values, df_y_pred.values))
# np.product(stack.shape)
# new_dim = (int(np.product(stack.shape)/2), 2)
# df_res = pd.DataFrame(np.reshape(stack, new_dim))

# dataset = "hiv_protease"
dataset = "ace_vaxinpad"

df_true = pd.read_csv(f"data/{dataset}/benchmark/ensemble/all_vs_all/group_1/y_true_cv_waac_aaindex_BUNA790103.csv", index_col=0)
df_pred = pd.read_csv(f"data/{dataset}/benchmark/ensemble/all_vs_all/group_1/y_pred_cv_waac_aaindex_BUNA790103.csv", index_col=0)

x, y = np.meshgrid(range(df_true.shape[1]),
                   range(df_true.shape[0]))

paths = glob(f"data/{dataset}/benchmark/ensemble/all_vs_all/group_1/y_pred_cv_*.csv")

z = pd.read_csv(paths[0], index_col=0).values
for p in paths[1:]:
    z += pd.read_csv(p, index_col=0).values

z = np.abs(df_true.values * len(paths) - z)
z /= len(paths)

# z = np.abs((df_true - df_pred).values)

source = pd.DataFrame({'x': x.ravel(),
                       'y': y.ravel(),
                       'z': z.ravel(),
                       "miss-classified": [">=50%" if i >= 0.5 else "<50%" for i in z.ravel()]})

alt.Chart(source).mark_rect().encode(
    x='x:O',
    y='y:O',
    color=alt.Color('miss-classified:N', scale=alt.Scale(domain=[">=50%", "<50%"], range=["black", "white"])),
    tooltip=["x", "y", "z"]
).properties(
    height=300,
    width=600
).save("chart.html")