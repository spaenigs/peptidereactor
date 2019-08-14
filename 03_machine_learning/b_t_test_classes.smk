import scripts.utils as utils
import pandas as pd
import numpy as np
import itertools


rule t_test_on_classification_error:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error.csv"
    run:
        from scipy import stats
        df = pd.read_csv(str(input), index_col=0)
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
        res.to_csv(str(output))


rule plot_t_test_on_classification_error:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error.csv"
    output:
        "00_data/out/neuropeptides/plots/{encoding}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error.pdf"
    script:
          "scripts/plot_t_test.py"


rule combine_cv_classes:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
                                 "{dataset}_{part}_normalized-{normalized}_cross_validation.csv",
                                 dataset=wildcards.dataset,
                                 part=wildcards.part,
                                 # TODO
                                 encoding=utils.PARAM_FREE_ENCODINGS + \
                                          utils.REST_ENCODINGS + utils.STRUC_ENCODINGS + \
                                          ["apaac", "paac",
                                           "geary", "moran", "nmbroto", "qsorder",
                                           "socnumber"],
                                 normalized=wildcards.normalized)
    output:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/classes/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    run:
        res = pd.DataFrame()
        for path in list(input):
            res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0, sort=True)
        res.to_csv(str(output))


rule all_vs_all_ttest_classes:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/classes/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/classes/" + \
         "{dataset}_{part}_normalized-{normalized}_ttest_all.csv"
    threads:
        8
    run:
        from pathos.multiprocessing import ProcessingPool as Pool
        from statsmodels.stats.multitest import multipletests
        from scipy import stats

        df = pd.read_csv(str(input), index_col=0)
        names = df.index
        res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                           index=names,
                           columns=names)

        def func(tuple):
            n1, n2 = tuple[0], tuple[1]
            errors_df1 = df.loc[n1, ~df.columns.isin(["train_size", "test_size"])]
            errors_df2 = df.loc[n2, ~df.columns.isin(["train_size", "test_size"])]
            m_diff = np.abs(np.mean(errors_df1 - errors_df2))
            sd = np.std(errors_df1 - errors_df2, ddof=1)
            N_training = df.loc[n1, "train_size"]
            N_testing = df.loc[n1, "test_size"]
            standard_error_mean_corr = sd * np.sqrt((1 / len(errors_df1)) + (N_testing / N_training))
            deg_freedom = len(errors_df1) - 1
            pv = 2 * stats.t.cdf(-m_diff / standard_error_mean_corr, deg_freedom)
            return n1, n2, 1.0 if np.isnan(pv) else pv

        p = Pool(threads)
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
        new_df.to_csv(str(output))


rule generate_network_data_from_classes:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/classes/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_all.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/classes/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_all.json"
    script:
         "scripts/generate_network_data.py"
