from scipy import stats
import pandas as pd
import numpy as np
import itertools

df = pd.read_csv(str(snakemake.input), index_col=0)
names = df.index
res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                   index=names,
                   columns=names)
for n1, n2 in itertools.combinations(names, 2):
    errors_df1 = df.loc[n1, ~df.columns.isin(["train_size", "test_size"])]
    errors_df2 = df.loc[n2, ~df.columns.isin(["train_size", "test_size"])]
    m_diff = np.abs(np.mean(errors_df1 - errors_df2))
    sd = np.std(errors_df1 - errors_df2, ddof=1)
    N_training = df.loc[n1, "train_size"]
    N_testing = df.loc[n1, "test_size"]
    standard_error_mean_corr = sd * np.sqrt((1 / len(errors_df1)) + (N_testing / N_training))
    deg_freedom = len(errors_df1) - 1
    res.loc[n1, n2] = 2 * stats.t.cdf(-m_diff / standard_error_mean_corr, deg_freedom)
res = res.transpose() * res
res.to_csv(str(snakemake.output))
