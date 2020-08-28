from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedKFold
from imblearn.under_sampling import RandomUnderSampler
from itertools import product

import pandas as pd
import numpy as np

from nodes.vis.single_dataset.scripts.utils import is_struc_based

dataset = "ace_vaxinpad"

df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)

indices = df_f1.apply(np.median)\
    .sort_values(ascending=False)\
    .index.tolist()

top_three_seq = indices[:3]
top_three_str = [i for i in indices if is_struc_based(i)][:3]

for encoding_pair in list(product(top_three_seq, top_three_str))[:1]:

    df1 = pd.read_csv(f"data/{dataset}/csv/all/{encoding_pair[0]}.csv", index_col=0)
    df2 = pd.read_csv(f"data/{dataset}/csv/all/{encoding_pair[1]}.csv", index_col=0)

    # justify training examples
    df1 = df1.loc[df2.index, :]

    X, y = df1.iloc[:, :-1].values, df1["y"]

    print(encoding_pair)

    ensemble_cv = {encoding_pair: {"ensemble_cv_scores": []}}
    for train_idx, test_idx in StratifiedKFold(n_splits=10).split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]

        res = {
            encoding_pair[0]: {"X": np.array([]), "y": np.array([]), "inner_cv_scores": []},
            encoding_pair[1]: {"X": np.array([]), "y": np.array([]), "inner_cv_scores": []}
        }

        for e in encoding_pair:

            df = pd.read_csv(f"data/{dataset}/csv/all/{e}.csv", index_col=0)
            X, y = df.iloc[:, :-1].values, df["y"]

            inner_cv_scores = []
            for train_idx_i, test_idx_i in StratifiedKFold(n_splits=10).split(X_train, y_train):
                X_train_i, y_train_i = X[train_idx_i], y[train_idx_i]
                X_test_i, y_test_i = X[test_idx_i], y[test_idx_i]

                # 0.3
                rus = RandomUnderSampler(sampling_strategy=1.0, random_state=42)
                X_train_i_resampled, y_train_i_resampled = rus.fit_resample(X_train_i, y_train_i)

                brf = RandomForestClassifier(n_jobs=12, random_state=42)
                brf.fit(X_train_i_resampled, y_train_i_resampled)
                y_pred_class = brf.predict(X_test_i)

                inner_cv_scores += [f1_score(y_true=y_test_i, y_pred=y_pred_class)]

            res[e]["inner_cv_scores"] = inner_cv_scores
            res[e]["X"], res[e]["y"] = X, y

        def get_proba(e):
            X, y = res[e]["X"], res[e]["y"]
            X_train_p, y_train_p = X[train_idx], y[train_idx]
            X_test_p, y_test_p = X[test_idx], y[test_idx]
            brf1 = RandomForestClassifier(n_jobs=12, random_state=42)
            brf1.fit(X_train_p, y_train_p)
            return brf1.predict_proba(X_test_p)[:, 1], X_test_p, y_test_p


        X1, y1 = df1.iloc[:, :-1].values, df1["y"]
        X1_test, y1_test = X1[test_idx], y1[test_idx]
        brf1 = RandomForestClassifier(n_jobs=12, random_state=42)
        brf1.fit(X1_test, y1_test)

        X2, y2 = df2.iloc[:, :-1].values, df2["y"]
        X2_test, y2_test = X2[test_idx], y2[test_idx]

        X, y = res[e]["X"][test_idx], res[e]["y"][test_idx]

        probas_model_1, y_test_model_1 = get_proba(encoding_pair[0])
        probas_model_2, _ = get_proba(encoding_pair[1])

        df_new = pd.DataFrame({"p1": probas_model_1, "p2": probas_model_2, "y": y_test_model_1})
        X, y = df_new.iloc[:, :-1].values, df_new["y"]

        X_train_e, y_train_e = X[train_idx], y[train_idx]
        X_test_e, y_test_e = X[test_idx], y[test_idx]

        ensemble = RandomForestClassifier(n_jobs=12, random_state=42)
        ensemble.fit(X_train_e, y_train_e)
        ypred = ensemble.predict(X_test_e)

        ensemble_cv[encoding_pair]["ensemble_cv_scores"] += [f1_score(y_test, ypred)]

    print(ensemble_cv)




