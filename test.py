import pandas as pd
import altair as alt
import numpy as np

from nodes.vis.single_dataset.scripts.utils import is_struc_based

ds = "ace_vaxinpad"
# ds = "hiv_protease"

df = pd.read_csv(f"data/{ds}/benchmark/metrics/f1.csv", index_col=0)
df_medians = df.apply(np.median).to_frame("median")

group = lambda enc: \
    "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:6]

df_medians["group"] = [group(x) for x in df_medians.index]
df_medians = df_medians.groupby("group").median()
df_medians["group"] = df_medians.index
df_medians.index = range(0, df_medians.shape[0])
df_medians["type"] = \
    ["structure based" if is_struc_based(e) else "sequence based" for e in df_medians["group"]]

df_time = pd.read_csv(f"data/{ds}/misc/benchmark/benchmark.csv", index_col=0)
df_time["meta_rule"] = df_time.index
df_time["group"] = [group(x) for x in df_time["meta_rule"]]

dur_tert_search = df_time.loc[df_time["meta_rule"]=="tertiary_structure_search", "s"][0]/60/60

# df_time.drop(df_time[df_time["category"] != "encoding"].index, inplace=True)

times = []
for n, s in df_medians.iterrows():
    sg = s["group"]
    df_tmp = df_time.loc[[g in sg and sg.startswith(g) for g in df_time["group"]], :]
    if sg == "aainde":
        df_tmp = df_time.loc[df_time["meta_rule"] == "encoding_aaindex", :]
    elif sg == "dist_f":
        df_tmp = df_time.loc[df_time["meta_rule"] == "distance_frequency", :]
    if not df_tmp.empty:
        times += [df_tmp["s"][0]]
    else:
        print(f"Could find {sg} in elapsed time df. Skipping.")
        df_medians.drop(df_medians.loc[df_medians["group"] == sg, :].index, inplace=True)

df_medians["time"] = [t/60/60 for t in times]
df_medians["time"] = df_medians\
    .apply(lambda s: s["time"] + dur_tert_search if is_struc_based(s.group) else s["time"], axis=1)\
    .to_list()

c = alt.Chart(df_medians).mark_point(filled=True).encode(
    x=alt.X("time:Q", scale=alt.Scale(type="log")),
    y="median:Q",
    color=alt.Color(
        "type:N",
        title="Encoding type",
        scale=alt.Scale(
            domain=["sequence based", "structure based"],
            range=["#7570b3", "#d95f02"]
        )
    ),
    tooltip=["group:N", "time:Q"]
)

text = c.mark_text(dy=10, size=9).encode(
    text="group:N"
)

(c + text).interactive().save("chart.html")