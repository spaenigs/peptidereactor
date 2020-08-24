import pandas as pd
import numpy as np
import altair as alt


def boxplot_altair(data, x, y, xtype='N', ytype='Q', size=40, width=400):


    """
    Python function to make boxplots in Altair
    """

    # Define variables and their types using f-strings in Python

    lower_box = f'q1({y}):{ytype}'
    # lower_whisker = f'min({y}):{ytype}'
    upper_box = f'q3({y}):{ytype}'
    # upper_whisker = f'max({y}):{ytype}'
    # median_whisker = f'median({y}):{ytype}'
    x_data = f'{x}:{xtype}'
    #
    # # lower plot
    # lower_plot = alt.Chart(data).mark_rule().encode(
    #     y=alt.Y(lower_whisker, axis=alt.Axis(title=y)),
    #     y2=lower_box,
    #     x=x_data
    # ).properties(
    #     width=width)
    #
    # # middle plot
    # middle_plot = alt.Chart(data).mark_bar(size=size).encode(
    #     y=lower_box,
    #     y2=upper_box,
    #     x=x_data
    # ).properties(
    #     width=width)
    #
    # # upper plot
    # upper_plot = alt.Chart(data).mark_rule().encode(
    #     y=upper_whisker,
    #     y2=upper_box,
    #     x=x_data
    # ).properties(
    #     width=width)
    #
    # # median marker line
    # middle_tick = alt.Chart(data).mark_tick(
    #     color='white',
    #     size=size
    # ).encode(
    #     y=median_whisker,
    #     x=x_data,
    # )
    #
    # # combine all the elements of boxplot to a single chart object
    # chart = lower_plot + middle_plot + upper_plot + middle_tick


    n = data["name"][0]
    data = data.loc[data["name"] == n, :]
    d_obj = data.loc[:, "fold_change"].describe()

    q1, median, q3 = d_obj["25%"], d_obj["50%"], d_obj["75%"]
    iqr = q3 - q1
    lower_outlier_cutoff, upper_outlier_cutoff = q1 - 1.5 * iqr, q3 + 1.5 * iqr

    data["q1"], data["median"], data["q3"] = q1, median, q3
    data["lower_outlier_cutoff"] = q1 - 1.5 * iqr
    data["upper_outlier_cutoff"] = q3 + 1.5 * iqr

    data["is_outlier"] = data[y].apply(lambda x: True if x < lower_outlier_cutoff else False)
    data["is_outlier"] = data[y].apply(lambda x: True if x > upper_outlier_cutoff else False)

    base = alt.Chart(data.loc[data["is_outlier"] == False, :])

    lower_rule = base.mark_rule().encode(
        x=x_data,
        y=f"min({y}):{ytype}",
        y2=f"q1({y}):{ytype}"
    ).properties(
        width=width
    )

    box = base.mark_bar(size=10).encode(
        x=f"{x}:{xtype}",
        y=f"q1({y}):{ytype}",
        y2=f"q3({y}):{ytype}"
    )

    median_tick = alt.Chart(data).mark_tick(
        color='white',
        size=size
    ).encode(
        y=f'min(median):Q',
        x=x_data,
    )

    upper_rule = base.mark_rule().encode(
        x=x_data,
        y=f"q3({y}):{ytype}",
        y2=f"max({y}):{ytype}"
    ).properties(
        width=width
    )

    base = alt.Chart(data.loc[data["is_outlier"], :])

    lower_outliers = base.mark_point().encode(
        x=f"{x}:{xtype}",
        y=f"{y}:{ytype}"
    )

    upper_outliers = base.mark_point().encode(
        x=f"{x}:{xtype}",
        y=f"{y}:{ytype}"
    )

    # return chart object
    return lower_rule + box + median_tick + upper_rule + upper_outliers + lower_outliers




df = pd.read_csv("intersections_200820.csv")
df["name"] = df["term_name"] + df["term_id"]

df["p_value"] = -np.log10(df["p_value"])
pvals = sorted(df["p_value"])

c = alt.Chart(df).mark_circle(stroke="black", strokeWidth=1.5).encode(
    x=alt.X("protein_enrichment:Q", title=None),
    y=alt.Y("name:N", sort="-x", title=None),
    size="protein_number:Q",
    color=alt.Color(
        "p_value:Q",
        title="-log10(p_value)",
        scale=alt.Scale(domain=[0, np.median(pvals), pvals[-1]], range=["black", "red", "white"], nice=20)
    ),
    row=alt.Row("source:N", title=None),
    tooltip="p_value:Q"
).properties(
    width=400,
    height=210
).resolve_scale(
    y="independent"
)

df_res = pd.DataFrame()
for n, s in df.iterrows():
    df_tmp = pd.DataFrame({"fold_change": s.filter(like="fold_change")})
    df_tmp.dropna(inplace=True)
    df_tmp["name"] = s["name"]
    df_tmp["source"] = s["source"]
    df_res = pd.concat([df_res, df_tmp])

bp = alt.Chart(df_res).mark_boxplot().encode(
    x=alt.X("fold_change:Q", title=None),
    y=alt.Y("name:N", title=None),
    row=alt.Row("source:N", title=None),
).properties(
    width=400,
    height=210
).resolve_scale(
    y="independent"
)


df_res = df_res.astype({'fold_change': 'float'})
bp2 = boxplot_altair(data=df_res, y="fold_change", ytype="Q", x="name", xtype="N", width=700, size=10)

alt.vconcat(
    alt.hconcat(
        c, bp
    ),
    bp2,
        config=alt.Config(
            axis=alt.AxisConfig(labelLimit=500)
        )
).save("test.html")

