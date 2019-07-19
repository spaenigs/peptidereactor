import pandas as pd
import numpy as np

df_tmp = pd.read_csv(str(snakemake.input), index_col=0, engine="c", usecols=lambda col: "y" not in col)

dfrows, dfcols = df_tmp.shape
dfzeros = np.count_nonzero(df_tmp == 0)
dftotal_eles = dfrows * dfcols

# TODO only filter by sparsity, AAINDEX: row/col ratio depends only on sequence
#  length and resulting interpolation
if (dfzeros / dftotal_eles) < snakemake.config["filter"]["max_sparse_ratio"]:
    df = pd.read_csv(str(snakemake.input), index_col=0, engine="c")
    df.to_csv(str(snakemake.output))
else:
    with open(str(snakemake.output), mode="w") as f:
        f.write("")
        f.flush()