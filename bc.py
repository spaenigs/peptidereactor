from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from imblearn.under_sampling import RandomUnderSampler
from itertools import product
from collections import Counter
from scipy import stats

import pandas as pd
import numpy as np

from nodes.vis.single_dataset.scripts.utils import is_struc_based

# dataset = "ace_vaxinpad"
dataset = "bce_bcm"
# dataset = "hiv_protease"

df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)

indices = df_f1.apply(np.median)\
    .sort_values(ascending=False)\
    .index.tolist()

top_three_seq = indices[:3]
top_three_str = [i for i in indices if is_struc_based(i)][:3]


def train_and_measure(X_train, y_train, X_test, y_test):
    brf = RandomForestClassifier(n_jobs=32, random_state=42)
    brf.fit(X_train, y_train)
    y_pred_class = brf.predict(X_test)
    return f1_score(y_true=y_test, y_pred=y_pred_class)


encoding_pairs = list(product(top_three_seq, top_three_str))[:1]
repeats = 1

df_total = pd.DataFrame()
for encoding_pair in encoding_pairs:

    df1 = pd.read_csv(f"data/{dataset}/csv/all/{encoding_pair[0]}.csv", index_col=0)
    df2 = pd.read_csv(f"data/{dataset}/csv/all/{encoding_pair[1]}.csv", index_col=0)

    # justify training examples
    df1 = df1.loc[df2.index, :]

    print(Counter(df1["y"]))

    X, y = df1.iloc[:, :-1].values, df1["y"]

    print(encoding_pair)

    for e in encoding_pair:

        df = pd.read_csv(f"data/{dataset}/csv/all/{e}.csv", index_col=0)
        X, y = df.iloc[:, :-1].values, df["y"]

        inner_cv_scores = []
        for train_idx, test_idx in RepeatedStratifiedKFold(n_repeats=repeats, n_splits=10).split(X, y):
            X_train, y_train = X[train_idx], y[train_idx]
            X_test, y_test = X[test_idx], y[test_idx]

            rus = RandomUnderSampler(sampling_strategy=0.3, random_state=42)
            X_train_i_resampled, y_train_i_resampled = rus.fit_resample(X_train, y_train)
            # print(Counter(y_train_i_resampled))
            f1s = train_and_measure(X_train_i_resampled, y_train_i_resampled, X_test, y_test)
            df_tmp = pd.DataFrame({"encoding": [f"{e}_resampled"], "f1": [f1s], "group": [str(encoding_pair)], "resampled": [True]})
            df_total = pd.concat([df_total, df_tmp])

            f1s = train_and_measure(X_train, y_train, X_test, y_test)
            df_tmp = pd.DataFrame({"encoding": [e], "f1": [f1s], "group": [str(encoding_pair)], "resampled": [False]})
            df_total = pd.concat([df_total, df_tmp])


def get_probas(df, train_idx, test_idx):
    X, y = df.iloc[:, :-1].values, df["y"]
    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]
    brf = RandomForestClassifier(n_jobs=32, random_state=42)
    brf.fit(X_train, y_train)
    return brf.predict_proba(X_test)[:, 1], y_test


def get_probas_resampled(df, train_idx, test_idx):
    X, y = df.iloc[:, :-1].values, df["y"]
    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]
    rus = RandomUnderSampler(sampling_strategy=0.3, random_state=42)
    X_train_i_resampled, y_train_i_resampled = rus.fit_resample(X_train, y_train)
    brf = RandomForestClassifier(n_jobs=32, random_state=42)
    brf.fit(X_train_i_resampled, y_train_i_resampled)
    return brf.predict_proba(X_test)[:, 1], y_test


# for encoding_pair in encoding_pairs:
#
#     print(encoding_pair)
#
#     df1 = pd.read_csv(f"data/{dataset}/csv/all/{encoding_pair[0]}.csv", index_col=0)
#     df2 = pd.read_csv(f"data/{dataset}/csv/all/{encoding_pair[1]}.csv", index_col=0)
#
#     # justify training examples
#     df1 = df1.loc[df2.index, :]
#
#     X, y = df1.iloc[:, :-1].values, df1["y"]
#
#     f1_scores = []
#     for train_idx, test_idx in RepeatedStratifiedKFold(n_repeats=repeats, n_splits=10).split(X, y):
#
#         X_train, y_train = X[train_idx], y[train_idx]
#         X_test, y_test = X[test_idx], y[test_idx]
#
#         df_res = pd.DataFrame()
#         for train_idx_i, test_idx_i in StratifiedKFold(n_splits=4).split(X_train, y_train):
#             probas1, y1_test = get_probas(df1, train_idx_i, test_idx_i)
#             probas2, y2_test = get_probas(df2, train_idx_i, test_idx_i)
#             df_new = pd.DataFrame({"p1": probas1, "p2": probas2, "y": y1_test})
#             df_res = pd.concat([df_res, df_new])
#
#         probasm1, _ = get_probas(df1, train_idx, test_idx)
#         probasm2, _ = get_probas(df2, train_idx, test_idx)
#
#         df_tmp = pd.DataFrame({"p1": probasm1, "p2": probasm2})
#
#         X_new, y_new = df_res.iloc[:, :-1].values, df_res["y"]
#
#         ensemble = RandomForestClassifier(n_jobs=32, random_state=42)
#         ensemble.fit(X_new, y_new)
#         ypred = ensemble.predict(df_tmp.values)
#
#         f1s = f1_score(y_test, ypred)
#
#         df_total = pd.concat([df_total, pd.DataFrame({
#             "encoding": [f"ensemble_{encoding_pair[0][:3]}_{encoding_pair[1][:3]}"],
#             "group": [str(encoding_pair)], "f1": [f1s], "resampled": [False]
#         })])
#
#         ### resampled
#
#         df_res = pd.DataFrame()
#         for train_idx_i, test_idx_i in StratifiedKFold(n_splits=4).split(X_train, y_train):
#             probas1, y1_test = get_probas_resampled(df1, train_idx_i, test_idx_i)
#             probas2, _ = get_probas_resampled(df2, train_idx_i, test_idx_i)
#             df_new = pd.DataFrame({"p1": probas1, "p2": probas2, "y": y1_test})
#             df_res = pd.concat([df_res, df_new])
#
#         probasm1, _ = get_probas(df1, train_idx, test_idx)
#         probasm2, _ = get_probas(df2, train_idx, test_idx)
#
#         df_tmp = pd.DataFrame({"p1": probasm1, "p2": probasm2})
#
#         X_new, y_new = df_res.iloc[:, :-1].values, df_res["y"]
#
#         ensemble = RandomForestClassifier(n_jobs=32, random_state=42)
#         ensemble.fit(X_new, y_new)
#         ypred = ensemble.predict(df_tmp.values)
#
#         f1s = f1_score(y_test, ypred)
#
#         df_total = pd.concat([df_total, pd.DataFrame({
#             "encoding": [f"ensemble_{encoding_pair[0][:3]}_{encoding_pair[1][:3]}_resampled"],
#             "group": [str(encoding_pair)], "f1": [f1s], "resampled": [True]
#         })])

df_total.to_csv(f"bc_data_{dataset}.csv")

