import pandas as pd
import altair as alt
import numpy as np


def scatter_div_chart(df_f1, df_div, dataset):

    def is_struc_based(e):
        if "asa" in e:
            return True
        elif "delaunay" in e:
            return True
        elif "disorder" in e:
            return True
        elif "hull" in e:
            return True
        elif "qsar" in e:
            return True
        elif "sse" in e:
            return True
        elif "ta" in e:
            return True
        else:
            return False

    def get_data(e1, e2):
        df1 = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/sequence_based/y_prob_cv_{e1}.csv",
                              index_col=0)
        df1_class = pd.read_csv(
            f"data/{dataset}/benchmark/ensemble/seq_vs_str/sequence_based/y_true_cv_{e1}.csv",
            index_col=0)
        df2 = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/structure_based/y_prob_cv_{e2}.csv",
                              index_col=0)
        return df1, df2, df1_class

    def ravel_and_annotate(df1, df2, df1_class, cat, div, e1, e2, f1_e1, f1_e2):
        df = pd.DataFrame({
            "x": df1.values.ravel(),
            "y": df2.values.ravel(),
            "class": df1_class.values.ravel()})
        df["diversity"] = \
            f"{cat} ({np.round(div, 2)}), e1: {e1} ({np.round(f1_e1, 1)}), e2: {e2} ({np.round(f1_e2, 1)})"
        return df


    medians = df_f1.apply(np.median).sort_values(ascending=False)
    indices = medians.index.tolist()
    indices_seq = indices[:15]
    indices_str = [i for i in indices if is_struc_based(i)][:5]

    df_div_sub = df_div.loc[indices_seq, indices_str]

    # crap
    i_crap, c_crap, div_crap = \
        "cgr_res_20_sf_0.5", "electrostatic_hull_6", df_div.loc["cgr_res_20_sf_0.5", "electrostatic_hull_6"]
    df1_crap, df2_crap, df1_crap_class = get_data(i_crap, c_crap)
    df_crap = ravel_and_annotate(df1_crap, df2_crap, df1_crap_class,
                                cat=f"(a) crap", div=div_crap, e1=i_crap, e2=c_crap,
                                f1_e1=medians[i_crap], f1_e2=medians[c_crap])

    # low
    i_low, c_low, div_low = sorted([(s, i, df_div.loc[s, i]) for i, s in df_div_sub.idxmin().items()], key=lambda x: x[2])[0]
    df1_low, df2_low, df1_low_class = get_data(i_low, c_low)
    df_low = ravel_and_annotate(df1_low, df2_low, df1_low_class,
                                cat=f"(b) low", div=div_low, e1=i_low, e2=c_low,
                                f1_e1=medians[i_low], f1_e2=medians[c_low])

    # mid
    i_mid, c_mid, div_mid = sorted([(i, j, df_div_sub.loc[i, j])
                                    for i, j in [(i, j) for i in df_div_sub.index for j in df_div_sub.columns]
                                    if 0.25 < df_div_sub.loc[i, j] < 0.26], key=lambda x: x[2])[0]
    df1_mid, df2_mid, df1_mid_class = get_data(i_mid, c_mid)
    df_mid = ravel_and_annotate(df1_mid, df2_mid, df1_mid_class,
                                cat=f"(c) mid", div=div_mid, e1=i_mid, e2=c_mid,
                                f1_e1=medians[i_mid], f1_e2=medians[c_mid])

    # high
    i_high, c_high, div_high = \
    sorted([(s, i, df_div.loc[s, i]) for i, s in df_div_sub.idxmax().items()], key=lambda x: x[2])[-1]
    df1_high, df2_high, df1_high_class = get_data(i_high, c_high)
    df_high = ravel_and_annotate(df1_high, df2_high, df1_high_class,
                                cat=f"(d) high", div=div_high, e1=i_high, e2=c_high,
                                f1_e1=medians[i_high], f1_e2=medians[c_high])

    # highest f1 each
    i_highest_f1, c_highest_f1 = indices_seq[0], indices_str[0]
    div_highest_f1 = df_div_sub.loc["ngram_a3_300", "delaunay_cartesian_product"]
    df1_highest_f1, df2_highest_f1, df1_highest_f1_class = get_data(i_highest_f1, c_highest_f1)
    df_highest_f1 = ravel_and_annotate(df1_highest_f1, df2_highest_f1, df1_highest_f1_class,
                                 cat=f"(e) highest f1", div=div_highest_f1, e1=i_highest_f1, e2=c_highest_f1,
                                 f1_e1=medians[i_highest_f1], f1_e2=medians[c_highest_f1])

    df = pd.concat([df_crap, df_low, df_mid, df_high, df_highest_f1])
    df = df.dropna()

    c = alt.Chart(df).mark_point(filled=True, size=5, opacity=1).encode(
        x=alt.X("x", title=f"Predicted probability e1"),
        y=alt.Y("y", title=f"Predicted probability e2"),
        color=alt.Color(
            "class:N",
            scale=alt.Scale(domain=[0, 1], range=["#7570b3", "#d95f02"])),
        column=alt.Column("diversity", header=alt.Header(labelOrient='left'))
    ).properties(
        width=250,
        height=250
    )

    return alt.hconcat(c, title=alt.TitleParams(text=["Increasing diversity w.r.t. performance", ""], anchor="middle"))



