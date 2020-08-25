import pandas as pd
import numpy as np
import altair as alt


def boxplot_altair(data, x, y, xtype, ytype, facet, sort, size=40, width=400):

    """
    Python function to make boxplots in Altair
    """

    # Define variables and their types using f-strings in Python

    # lower_box = f'q1({y}):{ytype}'
    # lower_whisker = f'min({y}):{ytype}'
    # upper_box = f'q3({y}):{ytype}'
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

    def describe_df(df):
        d_obj = df.loc[:, y].describe()
        q1, median, q3 = d_obj["25%"], d_obj["50%"], d_obj["75%"]
        iqr = q3 - q1
        lower_outlier_cutoff, upper_outlier_cutoff = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        df["q1"], df["median"], df["q3"] = q1, median, q3
        df["lower_outlier_cutoff"] = q1 - 1.5 * iqr
        df["upper_outlier_cutoff"] = q3 + 1.5 * iqr
        df["is_outlier"] = df[y].apply(lambda x: True if x < lower_outlier_cutoff else False)
        df["is_outlier"] = df[y].apply(lambda x: True if x > upper_outlier_cutoff else False)
        return df

    df_res = pd.DataFrame()
    for n in data[x].unique():
        df = data.loc[data[x]==n, :]
        df = describe_df(df)
        df_res = pd.concat([df_res, df])

    print(df_res.head())

    data = df_res

    field = sort

    charts = []
    for s in data[facet].unique():
        print(s)

        base = alt.Chart(data.loc[data["is_outlier"] == False, :]).transform_filter(
            alt.datum.source == s
        )

        lower_rule = base.mark_rule().encode(
            x=alt.X(x_data, sort=alt.EncodingSortField(field=field)),
            y=f"min({y}):{ytype}",
            y2=f"q1({y}):{ytype}"
        ).properties(
            width=width
        )

        box = base.mark_bar(size=10).encode(
            x=alt.X(x_data, sort=alt.EncodingSortField(field=field)),
            y=f"q1({y}):{ytype}",
            y2=f"q3({y}):{ytype}"
        )

        median_tick = alt.Chart(data).mark_tick(
            color='white',
            size=size
        ).encode(
            x=alt.X(x_data, sort=alt.EncodingSortField(field=field)),
            y=f'min(median):Q'
        ).transform_filter(
            alt.datum.source == s
        )

        upper_rule = base.mark_rule().encode(
            x=alt.X(x_data, sort=alt.EncodingSortField(field=field)),
            y=f"q3({y}):{ytype}",
            y2=f"max({y}):{ytype}"
        ).properties(
            width=width
        )

        base = alt.Chart(data.loc[data["is_outlier"], :]).transform_filter(
            alt.datum.source == s
        )

        lower_outliers = base.mark_point().encode(
            x=alt.X(x_data, sort=alt.EncodingSortField(field=field)),
            y=f"{y}:{ytype}"
        )

        upper_outliers = base.mark_point().encode(
            x=alt.X(x_data, sort=alt.EncodingSortField(field=field)),
            y=f"{y}:{ytype}"
        )

        c = alt.layer(lower_rule, box, median_tick, upper_rule, upper_outliers, lower_outliers)
        charts += [c]

    return alt.vconcat(*charts)


df = pd.DataFrame({
    "name": ["a"] * 10 + ["b"] * 10 + ["c"] * 10,
    "fold_change": list(np.random.random(9)) + [5] + list(np.random.random(20)),
    "source": ["foobar"] * 15 + ["tobos"] * 15
})

df = pd.read_json("data/temp/final/ace_vaxinpad/metrics/f1_box_plot_data.json")
print(df.head())

bp2 = boxplot_altair(data=df, y="Value", ytype="Q", x="Encoding", xtype="N", facet="metric", sort="median", width=700, size=10)

bp2.save("test.html")