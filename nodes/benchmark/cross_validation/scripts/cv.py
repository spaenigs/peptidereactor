from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold

import numpy as np
import pandas as pd


def get_splits(dataset_size):
    return int(np.round(dataset_size / (0.2 * dataset_size)))


def append_values(y, split_id):
    df = pd.DataFrame({f"split_{split_id}": y}).transpose()
    df.columns = [f"y_{i}" for i in range(df.shape[1])]
    return df


df = pd.read_csv(snakemake.input[0], index_col=0)
X, y = df.iloc[:, :-1].values, df["y"].values

brf = BalancedRandomForestClassifier(n_estimators=100, random_state=0)

cv = RepeatedStratifiedKFold(n_splits=get_splits(df.shape[0]), n_repeats=10)

df_y_true, df_y_pred, df_y_prob, df_imp = \
    pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

for i, (train_index, test_index) in enumerate(cv.split(X, y)):
    X_train, y_train = X[train_index], y[train_index]
    X_test, y_test = X[test_index], y[test_index]

    df_y_true = pd.concat([df_y_true, append_values(y_test, i)])

    brf.fit(X_train, y_train)

    df_imp_tmp = pd.DataFrame({f"fold_{i}": brf.feature_importances_})
    df_imp = pd.concat([df_imp, df_imp_tmp.transpose()])

    y_pred_class = brf.predict(X_test)
    df_y_pred = pd.concat([df_y_pred, append_values(y_pred_class, i)])

    y_pred_proba = brf.predict_proba(X_test)
    df_y_prob = pd.concat([df_y_prob, append_values(y_pred_proba[:, 1], i)])

df_y_true.to_csv(snakemake.output[0])
df_y_pred.to_csv(snakemake.output[1])
df_y_prob.to_csv(snakemake.output[2])
df_imp.to_csv(snakemake.output[3])
