import pandas as pd
import numpy as np
import altair as alt

df = pd.read_csv("intersections_200820.csv")
# df["term_name"] = df["term_name"].astype(str)
df2 = pd.read_csv("data.csv")
# df2["term_name"] = df2["term_name"].astype(str)
df = df.set_index('term_name').join(df2.set_index('term_name'), lsuffix="_")

df["term_name"] = df.index
df["name"] = df["term_name"] + " " + df["term_id"]

df["p_value"] = -np.log10(df["p_value"])
pvals = sorted(df["p_value"])

c = alt.Chart(df).mark_circle(stroke="black", strokeWidth=1.5).encode(
    x=alt.X("protein_enrichment:Q", title="Protein enrichment"),
    y=alt.Y("name:N", sort="-x", title=None),
    size=alt.Size("protein_number:Q", title="Protein number"),
    color=alt.Color(
        "p_value:Q",
        title="-log10(p-value)",
        scale=alt.Scale(domain=[0, np.median(pvals), pvals[-1]], range=["black", "red", "white"], nice=20)
    ),
    row=alt.Row("source_:N", title=None),
    tooltip="p_value:Q"
).properties(
    width=400,
    height=210
).resolve_scale(
    y="independent"
)

df_res = pd.DataFrame()
for n, s in df.iterrows():
    ad, ps = s["AD assoziiert"], s["beteiligte Ptotein_size"]
    df_tmp = pd.DataFrame({"fold_change": s.filter(like="fold_change")})
    df_tmp.dropna(inplace=True)
    df_tmp["name"] = s["name"] + f" ({ad}/{ps})"
    df_tmp["source_"] = s["source_"]
    df_res = pd.concat([df_res, df_tmp])

bp = alt.Chart(df_res).mark_boxplot().encode(
    x=alt.X("fold_change:Q", title="Fold enrichment"),
    y=alt.Y("name:N", title=None),
    row=alt.Row("source_:N", title=None),
).properties(
    width=400,
    height=210
).resolve_scale(
    y="independent"
)

t = alt.Chart(df_res).mark_text().encode(
    text="ratio:N"
)

alt.hconcat(
    c, bp,
    config=alt.Config(
        axis=alt.AxisConfig(labelLimit=500)
    )
).resolve_scale(
    color="independent",
    size="independent"
).save("test.html")

