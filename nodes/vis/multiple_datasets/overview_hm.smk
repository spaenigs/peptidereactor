from functools import reduce

import pandas as pd
import altair as alt
import numpy as np

import re

TOKEN = config["token"]


def compute_median_f1(path):
    dataset = re.findall("data/(.*?)/", path)[0]
    return pd\
        .read_csv(path, index_col=0)\
        .apply(np.median)\
        .to_frame(dataset)


rule overview_hm:
    input:
        config["benchmark_dirs_in"]
    output:
        config["html_dir_out"] + f"data/overview_hm/data.json",
        temp(f"data/temp/{TOKEN}/overview_hm.json")
    run:
        paths = [p + "metrics/f1.csv" for p in list(input)]

        df_res = pd.concat(
            axis=1,
            objs=reduce(
                lambda dfs, p: dfs + [compute_median_f1(p)],
                paths, [])).\
            sort_index().\
            fillna(0.0).\
            T\
            .iloc[:, :100]

        x, y = np.meshgrid(df_res.columns, df_res.index)

        source = pd.DataFrame({"Encoding": x.ravel(),
                               "Dataset": y.ravel(),
                               "F1": df_res.values.ravel()
                               })

        l = np.hstack((np.array([df_res.index]).T, df_res.values))
        l = sorted(l, key=lambda x: [*x[1:]], reverse=True)
        sorted_datasets = np.array(l)[:, 0]

        url = output[0].replace(config["html_dir_out"], "")

        chart_json = alt.Chart(url).mark_rect().encode(
            x='Encoding:N',
            y=alt.Y(
                'Dataset:N',
                sort=alt.Sort(alt.SortArray(sorted_datasets))
            ),
            color='F1:Q',
            tooltip=["Encoding:N", "Dataset:N", "F1:Q"]
        ).properties(
            width=900
        ).to_json(
            indent=None
        )

        source.to_json(output[0], orient="records")

        with open(output[1], "w") as f:
             f.write(chart_json)
             f.flush()