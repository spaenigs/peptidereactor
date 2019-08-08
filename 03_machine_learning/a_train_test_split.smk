from pandas.errors import EmptyDataError
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import f1_score

import scripts.utils as utils
from sklearn.model_selection import train_test_split, KFold
import pandas as pd
import numpy as np
import itertools

rule split_train_test:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
         "{dataset}_{part}_{type}.csv"
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/splitted/train/" + \
          "{dataset}_{part}_{type}_normalized-{normalized}.csv",
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/splitted/test/" + \
          "{dataset}_{part}_{type}_normalized-{normalized}.csv"
    run:
        df = pd.read_csv(str(input), index_col=0)
        df_train, df_test = train_test_split(df, test_size=0.1, random_state=42)
        df_train.to_csv(str(output[0]))
        df_test.to_csv(str(output[1]))


rule cross_validation:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/splitted/train/" + \
         "{dataset}_{part}_{type}_normalized-{normalized}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
        "{dataset}_{part}_{type}_normalized-{normalized}.csv"
    run:
        df = pd.read_csv(str(input), index_col=0)
        X, y = df.iloc[:, :-1].values, df["y"]
        res = {wildcards.type: {}}
        cnt = 1
        for train_index, test_index in KFold(n_splits=10, random_state=42).split(X):
            X_train, X_validation = X[train_index], X[test_index]
            y_train, y_validation = y[train_index], y[test_index]
            clfLDA = LinearDiscriminantAnalysis()
            clfLDA.fit(X_train, y_train)
            predsLDA = clfLDA.predict(X_validation)
            res[wildcards.type]["run_" + str(cnt)] = np.around((1 - f1_score(y_validation, predsLDA)) * 100, decimals=1)
            cnt += 1
        res[wildcards.type]["train_size"] = X_train.shape[0]
        res[wildcards.type]["test_size"] = X_validation.shape[0]
        pd.DataFrame(res).to_csv(str(output))


# def get_input(wildcards):
#     if wildcards.encoding in \
#             utils.PARAM_BASED_ENCODINGS + \
#             utils.REST_ENCODINGS:  # aaindex, psekraac
#         return expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
#                       "{dataset}_{part}_{type}_normalized-{normalized}.csv",
#                       dataset=wildcards.dataset,
#                       part=wildcards.part,
#                       encoding=wildcards.encoding,
#                       type=utils.get_type(wildcards.encoding, config),
#                       normalized=wildcards.normalized)
#     else:
#         return expand("00_data/out/neuropeptides/neuropeptides_ds1/encodings/{encoding}/part/" + \
#                       "csv/final/geom_median/tsne/normalized-{normalized}/{dataset}_{part}_{type}.csv",
#                       dataset=wildcards.dataset,
#                       part=wildcards.part,
#                       encoding=utils.PARAM_FREE_ENCODINGS + utils.STRUC_ENCODINGS,
#                       type=utils.get_type(wildcards.encoding, config),
#                       normalized=wildcards.normalized)

rule collect_cross_validation:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
                                 "{dataset}_{part}_{type}_normalized-{normalized}.csv",
                                 dataset=wildcards.dataset,
                                 part=wildcards.part,
                                 encoding=wildcards.encoding,
                                 type=utils.get_type(wildcards.encoding, config),
                                 normalized=wildcards.normalized)
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    run:
        df_res = pd.DataFrame()
        for path in list(input):
            try:
                df = pd.read_csv(path, index_col=0)
                df_res = pd.concat([df_res, df], axis=1, sort=True)
            except EmptyDataError:
                pass
        df_res = df_res.transpose()
        df_res.to_csv(str(output))  # reindex(sorted(df_res.columns, key=lambda x: int(x[4:])), axis=1)


rule t_test_on_classification_error_single:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error.csv"
    script:
        "scripts/t_test_on_classification_error.py"


rule combine_encodings_for_t_test:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
                                 "{dataset}_{part}_{type}_normalized-{normalized}.csv",
                                 dataset=wildcards.dataset,
                                 part=wildcards.part,
                                 encoding=utils.PARAM_FREE_ENCODINGS + utils.STRUC_ENCODINGS,
                                 type=utils.get_type(wildcards.encoding, config),
                                 normalized=wildcards.normalized)
    output:
        "00_data/out/{dataset}/{dataset}_{part}/cv_combined/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    run:
        res = pd.DataFrame()
        for path in list(input):
            res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0)
        res.to_csv(str(output))

"""
00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/{dataset}_{part}_normalized-{normalized}_cross_validation.csv
00_data/out/{dataset}/{dataset}_{part}/            cv_combined/{dataset}_{part}_normalized-{normalized}_cross_validation.csv
"""
rule t_test_on_classification_error_aggregated:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/cv_combined/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/cv_combined/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error.csv"
    script:
        "scripts/t_test_on_classification_error.py"


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