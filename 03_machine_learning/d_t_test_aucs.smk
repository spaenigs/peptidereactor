import scripts.utils as utils
import pandas as pd
import numpy as np
import itertools


rule t_test_on_aucs:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
         "{dataset}_{part}_normalized-{normalized}_cross_validation_aucs.csv"
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
          "{dataset}_{part}_normalized-{normalized}_ttest_aucs.csv"
    run:
        from scipy import stats
        df = pd.read_csv(str(input), index_col=0)
        names = df.index
        res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                           index=names,
                           columns=names)
        for n1, n2 in itertools.combinations(names, 2):
            _, pvalue = stats.ttest_rel(df.loc[n1, :].values, df.loc[n2, :].values)
            res.loc[n1, n2] = pvalue
        res = res.transpose() * res
        res.to_csv(str(output))


rule plot_t_test_on_aucs:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_aucs.csv"
    output:
        "00_data/out/neuropeptides/plots/{encoding}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_aucs.pdf"
    script:
          "scripts/plot_t_test.py"


rule combine_cv_aucs:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
                                 "{dataset}_{part}_normalized-{normalized}_cross_validation_aucs.csv",
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
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/aucs/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    run:
        res = pd.DataFrame()
        for path in list(input):
            res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0, sort=True)
        res.to_csv(str(output))


rule all_vs_all_ttest_aucs:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/aucs/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/aucs/" + \
         "{dataset}_{part}_normalized-{normalized}_ttest_all.csv"
    threads:
        8
    run:
        from pathos.multiprocessing import ProcessingPool as Pool
        from scipy import stats

        df = pd.read_csv(str(input), index_col=0)
        names = df.index
        res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                           index=names,
                           columns=names)

        def func(tuple):
            n1, n2 = tuple[0], tuple[1]
            _, pv = stats.ttest_rel(df.loc[n1, :].values, df.loc[n2, :].values)
            return n1, n2, 1.0 if np.isnan(pv) else pv

        p = Pool(threads)
        pvs = p.map(func, list(itertools.combinations(res.index, 2)))

        for n1, n2, pv in pvs:
            res.loc[n1, n2] = pv

        res = res.transpose() * res

        res.to_csv(str(output))


rule generate_network_data_from_aucs:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/aucs/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_all.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/aucs/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_all.json"
    script:
         "scripts/generate_network_data.py"
