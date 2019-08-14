from typing import List

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


def sort_df_to_lower_tril(df):
    """
    1) Sort rows
        a   b    c              a   b  c
    a  NaN NaN  1.0         b  NaN NaN  NaN
    b  NaN NaN  NaN  ==>    a  NaN NaN  1.0
    c  3.0 NaN  2.0         c  3.0 NaN  2.0
    2) Sort cols
        a   b  c                c   a   b
    b  NaN NaN  NaN         b  NaN NaN NaN
    a  NaN NaN  1.0  ==>    a  1.0 NaN NaN
    c  3.0 NaN  2.0         c  2.0 3.0 NaN

    """
    sorted_row_names = np.array(sorted([[row[0], ~np.isnan(row[1]).sum()]
                                        for row in df.iterrows()],
                                       key=lambda tup: tup[1]))[:, 0]
    df = df.reindex(sorted_row_names)

    sorted_col_names = np.array(sorted([[row[0], np.isnan(row[1]).sum()]
                                        for row in df.transpose().iterrows()],
                                       key=lambda tup: tup[1]))[:, 0]
    return df[sorted_col_names]


def concat(paths: List[str]) -> pd.DataFrame:
    df_res = pd.DataFrame()
    for path in paths:
        try:
            df = pd.read_csv(path, index_col=0)
            df_res = pd.concat([df_res, df], axis=1, sort=True)
        except EmptyDataError:
            pass
    return df_res if df_res.empty else sort_df_to_lower_tril(df_res)


concat(list(snakemake.input)).to_csv(str(snakemake.output))