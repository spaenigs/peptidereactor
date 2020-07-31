import pandas as pd
import numpy as np
import altair as alt

df = pd.read_csv("data/ace_vaxinpad/benchmark/metrics/f1.csv", index_col=0)
df_medians = df.apply(np.median).to_frame("median")
group = lambda enc: "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:6]
df_medians["group"] = [group(x) for x in df_medians.index]

alt.layer(
    alt.Chart(df_medians).mark_circle().encode(
        x=alt.X('group:N', axis=alt.Axis(grid=True)),
        y=alt.Y('median:Q', bin=True),
        size='count()',
        tooltip="count()"
    )
).save("chart.html")