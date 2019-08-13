from pandas.errors import EmptyDataError
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split, KFold

import scripts.utils as utils
import pandas as pd
import numpy as np


def get_input_split_train_test(wildcards):
    if wildcards.encoding in utils.REST_ENCODINGS:
        return "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
               "geom_median/tsne/normalized-yes/{dataset}_{part}_{type}.csv"
    else:
        return "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
               "{dataset}_{part}_{type}.csv"

rule split_train_test:
    input:
         get_input_split_train_test
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
        "{dataset}_{part}_{type}_normalized-{normalized}.csv",
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
        "{dataset}_{part}_{type}_normalized-{normalized}_scores.csv"
    run:
        df = pd.read_csv(str(input), index_col=0)
        X, y = df.iloc[:, :-1].values, df["y"]

        res, scores = {wildcards.type: {}}, {wildcards.type: np.array([])}

        cnt = 1
        for train_index, test_index in KFold(n_splits=10, random_state=42).split(X):
            X_train, X_validation = X[train_index], X[test_index]
            y_train, y_validation = y[train_index], y[test_index]
            clfLDA = LinearDiscriminantAnalysis()
            clfLDA.fit(X_train, y_train)

            predsLDA_classes = clfLDA.predict(X_validation)
            res[wildcards.type]["run_" + str(cnt)] = \
                np.around((1 - f1_score(y_validation, predsLDA_classes)) * 100, decimals=1)

            predsLDA_scores = clfLDA.predict_proba(X_validation)
            scores[wildcards.type] = \
                np.append(scores[wildcards.type], np.round(1 - predsLDA_scores[:, 1], 2))

            cnt += 1

        res[wildcards.type]["train_size"] = X_train.shape[0]
        res[wildcards.type]["test_size"] = X_validation.shape[0]
        pd.DataFrame(res).to_csv(str(output[0]))

        pd.DataFrame(scores).to_csv(str(output[1]))


def get_single_type(wildcards):
    type_, = glob_wildcards(f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/" + \
                            f"encodings/{wildcards.encoding}/" + \
                            f"csv/final/geom_median/tsne/normalized-{wildcards.normalized}/" + \
                            f"{wildcards.dataset}_{wildcards.part}_{{type}}.csv")
    return type_

rule collect_cross_validation_classes:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
                                 "{dataset}_{part}_{type}_normalized-{normalized}.csv",
                                 dataset=wildcards.dataset,
                                 part=wildcards.part,
                                 encoding=wildcards.encoding,
                                 type=get_single_type(wildcards) \
                                     if wildcards.encoding in utils.REST_ENCODINGS \
                                     else utils.get_type(wildcards.encoding, config),
                                 normalized=wildcards.normalized),
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    run:
        def concat_dfs(list_in):
            df_res = pd.DataFrame()
            for path in list_in:
                try:
                    df = pd.read_csv(path, index_col=0)
                    df_res = pd.concat([df_res, df], axis=1, sort=True)
                except EmptyDataError:
                    pass
            return df_res.transpose()
        concat_dfs(list(input)).to_csv(str(output))


rule collect_cross_validation_scores:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/part/" + \
                                 "{dataset}_{part}_{type}_normalized-{normalized}_scores.csv",
                                 dataset=wildcards.dataset,
                                 part=wildcards.part,
                                 encoding=wildcards.encoding,
                                 type=get_single_type(wildcards) \
                                     if wildcards.encoding in utils.REST_ENCODINGS \
                                     else utils.get_type(wildcards.encoding, config),
                                 normalized=wildcards.normalized)
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_scores.csv"
    run:
        def concat_dfs(list_in):
            df_res = pd.DataFrame()
            for path in list_in:
                try:
                    df = pd.read_csv(path, index_col=0)
                    df_res = pd.concat([df_res, df], axis=1, sort=True)
                except EmptyDataError:
                    pass
            return df_res.transpose()
        concat_dfs(list(input)).to_csv(str(output))
