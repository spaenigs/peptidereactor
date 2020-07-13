import altair as alt


def elapsed_time_chart(df_time):
    df_time["meta_rule"] = df_time.index
    c = alt.Chart(df_time).mark_point(opacity=1.0, size=50).encode(
        x=alt.X("meta_rule:N", sort="-y", title="Meta rule"),
        y=alt.Y("time_in_hours:Q", scale=alt.Scale(type='log'),  title="Hours (log scale)"),
        color=alt.Color("class:N", scale=alt.Scale(domain=["encoding", "benchmark", "utils"], range=["#1b9e77", "#d95f02", "#7570b3"])),
        shape=alt.Shape("class:N", scale=alt.Scale(domain=["encoding", "benchmark", "utils"], range=["cross", "square", "triangle-right"])),
        tooltip=["meta_rule", "time_in_hours:Q", "class:N"]
    ).transform_calculate(
        time_in_hours="datum.s/60/60"
    ).add_selection(
        alt.selection_single()
    ).properties(width=900)
    return alt.hconcat(c, title=alt.TitleParams(text=["Elapsed time", ""], anchor="middle"))