import glob
import itertools
import os
import re

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
from scipy import interpolate, stats


def write_empty_file(path):
    with open(path, mode="w") as f:
        f.write("")
        f.flush()

def interpolate_to(X: np.array, dim: int):
    ydim = X.shape[1]
    x = np.arange(0, ydim if ydim >= 2 else ydim + 1)
    func = interpolate.interp1d(x, X, kind="linear")
    return func(np.linspace(0, ydim - 1, dim))


dataset = snakemake.wildcards.dataset
part = snakemake.wildcards.part
normalized = snakemake.wildcards.normalized
encoding = snakemake.wildcards.encoding

try:

    actual_csv = str(snakemake.input[0])
    ds1 = pd.read_csv(actual_csv, engine="c", index_col=0).iloc[:, :-1].astype("float64").values

    files = [f for f in
             glob.glob(
                 f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" +
                 f"{dataset}_{part}_ifeature_*_normalized-{normalized}*")
             if os.path.getsize(f) > 0]

    files_filtered = \
        list(map(lambda tup: tup[1], filter(lambda tup: tup[0] == actual_csv, itertools.combinations(files, 2))))

    pattern = f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
              f"{dataset}_{part}_ifeature_(.*?)_normalized-{normalized}.csv"

    row_name = re.match(pattern, actual_csv).group(1)

    inner_dict = {}
    for _f in files:
        inner_dict.setdefault(re.match(pattern, _f).group(1), "NA")

    res_corr = {row_name: inner_dict}

    if len(files_filtered) is not 0:
        for ff in files_filtered:
            col_name = re.match(pattern, ff).group(1)
            ds2 = pd.read_csv(ff, engine="c", index_col=0).iloc[:, :-1].values
            to_dim = int(np.ceil(np.mean([ds1.shape[1], ds2.shape[1]])))
            ds1_interpolated, ds2_interpolated = interpolate_to(ds1, to_dim), interpolate_to(ds2, to_dim)
            corr, _ = stats.pearsonr(ds1_interpolated.ravel(), ds2_interpolated.ravel())
            res_corr[row_name][col_name] = 1 - corr

    df = pd.DataFrame(res_corr)
    df.to_csv(str(snakemake.output))


except EmptyDataError as e:
            write_empty_file(str(output))