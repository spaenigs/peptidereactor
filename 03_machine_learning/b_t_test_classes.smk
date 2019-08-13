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
    run:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        res = pd.read_csv(str(input), index_col=0)

        # https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
        fig, ax = plt.subplots()
        pc = ax.imshow(res.values, vmin=0, vmax=1, cmap='Greys')
        ax.set_xticks(np.arange(len(res.columns)))
        ax.set_yticks(np.arange(len(res.index)))
        ax.set_xticklabels(res.columns)
        ax.set_yticklabels(res.index)

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        cbar = ax.figure.colorbar(pc)
        cbar.ax.set_ylabel("p-value, (*): <= 0.05, (**): <= 0.01", rotation=-90, va="bottom")

        for i in range(len(res.index)):
            for j in range(len(res.columns)):
                if res.loc[res.index[i], res.columns[j]] <= 0.01:
                    ax.text(j, i, "(**)", ha="center", va="center", color="black")
                elif res.loc[res.index[i], res.columns[j]] <= 0.05:
                    ax.text(j, i, "(*)", ha="center", va="center", color="black")

        fig.tight_layout()
        fig.set_size_inches(9, 7)

        plt.title(f"{utils.get_encoding_description(wildcards.encoding)} ({wildcards.encoding.upper()}):\n" +
                  "paired t-test with corrected variance on classification error")

        plt.savefig(str(output), bbox_inches="tight")


rule combine:
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
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    run:
        res = pd.DataFrame()
        for path in list(input):
            res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0, sort=True)
        res.to_csv(str(output))


rule all_vs_all_ttest:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/" + \
         "{dataset}_{part}_normalized-{normalized}_ttest_error_all.csv"
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
        pvs = p.map(func, list(itertools.combinations(res.index, 2)))

        for n1, n2, pv in pvs:
            res.loc[n1, n2] = pv

        res = res.transpose() * res

        res.to_csv(str(output))


rule generate_network_data:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error_all.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error_all.json"
    run:
        import json
        import re
        df = pd.read_csv(str(input), index_col=0)

        # TODO blosum62
        nodes = [{"id": idx, "group": re.match("(.*?)(\d{1,2}[ABC]?)?encoder.*", idx).group(1)}
                 for idx in df.index]

        groups = set(list(map(lambda idx: re.match("(.*?)(type\d{1,2}[ABC]?)?encoder.*", idx).group(1), df.index)))

        links = []
        for g1, g2 in itertools.combinations(groups, 2):
            df_filtered = df.filter(regex=f"^{g1}(type\d{{1,2}}[ABC]?)?encoder.*", axis=0)  # filter by rows
            df_filtered = df_filtered.filter(regex=f"^{g2}(type\d{{1,2}}[ABC]?)?encoder.*")  # filter by columns
            min_index_idx, min_col_idx, val = \
                sorted([(n1, n2, df_filtered.loc[n1, n2])
                        for n1, n2 in zip(df_filtered.index, df_filtered.idxmin(axis=1).values)],
                       key=lambda triple: triple[2])[0]
            if val <= 0.05:
                links += [{"source": min_index_idx, "target": min_col_idx, "pvalue": val, "within": 0}]

        for group in groups:
            df_f2 = df.filter(regex=f"^{group}(type\d{{1,2}}[ABC]?)?encoder.*", axis=0)
            df_f2 = df_f2.filter(regex=f"^{group}(type\d{{1,2}}[ABC]?)?encoder.*", axis=1)
            if np.sum(df_f2.shape) > 2:
                for n1, n2 in itertools.combinations(df_f2.index, 2):
                    links += [{"source": n1, "target": n2, "pvalue": df_f2.loc[n1, n2], "within": 1}]

        res = json.dumps({"nodes": nodes, "links": links})

        with open(str(output), mode="a") as f:
            f.write(res)
            f.flush()