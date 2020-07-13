import altair as alt
import pandas as pd
import numpy as np


def par_coord_chart(df_f1, df_fir):

    df_fir.index = [e.replace(".txt", "") for e in df_fir.index]
    df_fir["Encoding"] = df_fir.index

    df_f1_top = df_f1.apply(np.median).sort_values(ascending=False).to_frame("Value").iloc[:142, :]
    df_f1_top["Type"] = "top"
    df_f1_flop = df_f1.apply(np.median).sort_values().to_frame("Value").iloc[:142, :]
    df_f1_flop["Type"] = "flop"
    df_f1_top = pd.concat([df_f1_top, df_f1_flop])
    df_f1_top["Encoding"] = df_f1_top.index
    df_f1_top["x"] = "F1"

    df_fir_top = df_fir.loc[df_f1_top.index, :]
    df_fir_top.columns = ["Value", "Encoding"]
    df_fir_top["Type"] = df_f1_top["Type"]
    df_fir_top["x"] = "FIR"

    source = pd.concat([df_f1_top, df_fir_top])

    d, r = ["top", "flop"], ["#7570b3", "#d95f02"]

    selection = alt.selection_single(on="mouseover")

    c = alt.Chart(source).mark_line().encode(
        x=alt.X("x:N", title=None),
        y="Value",
        detail="Encoding:N",
        color=alt.condition(
            selection,
            alt.Color("Type:N", scale=alt.Scale(domain=d, range=r)),
            alt.value("lightgrey")),
        tooltip=["Encoding"]
    ).add_selection(
        selection
    ).properties(width=900)

    return alt.hconcat(c, title=alt.TitleParams(text=["Case study: feature importance ratio", ""], anchor="middle"))
