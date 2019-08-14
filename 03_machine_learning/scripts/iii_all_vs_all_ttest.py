from pathos.multiprocessing import ProcessingPool as Pool
from statsmodels.stats.multitest import multipletests
from scipy import stats

import pandas as pd
import numpy as np
import itertools

df = pd.read_csv(str(snakemake.input), index_col=0)
names = df.index
res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                   index=names,
                   columns=names)


def func(tuple):
    n1, n2 = tuple[0], tuple[1]
    if snakemake.wildcards.computation == "classes":
        errors_df1 = df.loc[n1, ~df.columns.isin(["train_size", "test_size"])]
        errors_df2 = df.loc[n2, ~df.columns.isin(["train_size", "test_size"])]
        m_diff = np.abs(np.mean(errors_df1 - errors_df2))
        sd = np.std(errors_df1 - errors_df2, ddof=1)
        N_training = df.loc[n1, "train_size"]
        N_testing = df.loc[n1, "test_size"]
        standard_error_mean_corr = sd * np.sqrt((1 / len(errors_df1)) + (N_testing / N_training))
        deg_freedom = len(errors_df1) - 1
        pv = 2 * stats.t.cdf(-m_diff / standard_error_mean_corr, deg_freedom)
        # return n1, n2, 1.0 if np.isnan(pv) else pv
    elif snakemake.wildcards.computation == "scores" or snakemake.wildcards.computation == "aucs":
        _, pv = stats.ttest_rel(df.loc[n1, :].values, df.loc[n2, :].values)
        # return n1, n2, 1.0 if np.isnan(pv) else pv
    else:
        raise ValueError(f"Unknown wildcard for computation: {snakemake.wildcards.computation}")
    return n1, n2, 1.0 if np.isnan(pv) else pv


p = Pool(snakemake.threads)
pvs_triple = p.map(func, list(itertools.combinations(res.index, 2)))

for n1, n2, pv in pvs_triple:
    res.loc[n1, n2] = pv

res = res.transpose() * res

idxs = np.tril_indices(len(res.index), k=-1)
pvs = res.values[idxs]

# p-value correction for multiple tests using Benjamini/Hochberg (Biostatistics p.90)
# if rejected = true, accept H1
reject, pvalues_corrected, _, _ = multipletests(pvs, method="fdr_by")

print(f"before: {sum(list(map(lambda t: t[2] <= 0.05, pvs_triple)))}")
print(f"after: {sum(reject)}")

new_df_array = np.ones(res.shape)
new_df_array[idxs] = pvalues_corrected

new_df = pd.DataFrame(new_df_array.transpose() * new_df_array)
new_df.index, new_df.columns = res.index, res.columns
new_df.to_csv(str(snakemake.output))
