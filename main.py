import pandas as pd
import numpy as np

from scipy import stats

import itertools
import re
import glob

from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.metrics import f1_score

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression

df1 = pd.read_csv("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/" +
                  "neuropeptides_ds1_apaacencoder_lambda-9.csv", index_col=0)

df2 = pd.read_csv("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/" +
                 "neuropeptides_ds1_apaacencoder_lambda-16.csv", index_col=0)


df3 = pd.read_csv("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/" +
                  "neuropeptides_ds1_apaacencoder_lambda-30.csv", index_col=0)

# X_global = df1.iloc[:, :-1].values
# y_global = df1["y"]
#
# X, X_test, y, y_test = train_test_split(X_global, y_global, test_size=0.1, random_state=42)
#
# print(f"Train: {X.shape} {y.shape}\n Test: {X_test.shape}  {y_test.shape}")
#
# kf = KFold(n_splits=10, random_state=2652124)

def get_classification_error(df):
    X_global = df.iloc[:, :-1].values
    y_global = df["y"]
    X, X_test, y, y_test = train_test_split(X_global, y_global, test_size=0.1, random_state=42)
    res = []
    train_size = []
    test_size = []
    for train_index, test_index in KFold(n_splits=10, random_state=42).split(X):
        X_train, X_validation = X[train_index], X[test_index]
        y_train, y_validation = y[train_index], y[test_index]
        clfLDA = LinearDiscriminantAnalysis()
        clfLDA.fit(X_train, y_train)
        predsLDA = clfLDA.predict(X_validation)
        res.append(np.around((1 - f1_score(y_validation, predsLDA)) * 100, decimals=1))
        train_size.append(X_train.shape[0])
        test_size.append(X_validation.shape[0])
    return np.array(res), train_size[0], test_size[0]

# print(f"E1: {get_classification_error(df1)}")
# print(f"E2: {get_classification_error(df2)}")
# print(f"E3: {get_classification_error(df3)}")

# Index= ['apaac_lambda-9', 'apaac_lambda-16', 'apaac_lambda-30']
# Cols = ['apaac_lambda-9', 'apaac_lambda-16', 'apaac_lambda-30']



def get_dataframes(paths):
    for path in paths:
        yield pd.read_csv(path, index_col=0),  path # ".*?(\d{1,2}).csv"# re.match(".*_(.*?encoder.*?).csv", path).group(1)

file_paths = sorted(glob.glob("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/*.csv"),
                    key=lambda x: int(re.match(".*?(\d{1,2}).csv", x).group(1)))

res = pd.DataFrame(np.zeros((30, 30)) + 1.0, index=file_paths, columns=file_paths)

for (df_tmp1, n1), (df_tmp2, n2) in itertools.combinations(get_dataframes(file_paths), 2):   # itertools.combinations([(df1, 'apaac_lambda-9'), (df2, 'apaac_lambda-16'), (df3, 'apaac_lambda-30')], 2):
    errors_e1, N_training, N_testing = get_classification_error(df_tmp1)  # np.array([7.4, 18.1, 13.7, 17.5, 13.0, 12.5, 8.9, 12.1, 12.4, 7.4]), 9, 1  #
    errors_e2, _, _ = get_classification_error(df_tmp2)  # np.array([9.9, 11.0, 5.7, 12.5, 2.7, 6.6, 10.6, 6.4, 12.5, 7.8]), 9, 1  #
    m_diff = np.abs(np.mean(errors_e1 - errors_e2))
    sd = np.std(errors_e1 - errors_e2, ddof=1)
    standard_error_mean_corr = sd * np.sqrt((1 / len(errors_e1)) + (N_testing / N_training))
    df = len(errors_e1) - 1
    res.loc[n1, n2] = 2 * stats.t.cdf(-m_diff / standard_error_mean_corr, df)
    # if res.loc[n1, n2] <= 0.05:
    #     print(res.loc[n1, n2])

import matplotlib.pyplot as plt
plt.pcolor(res.values * res.values.transpose(), vmin=0, vmax=1, cmap='Purples')
plt.xticks(np.arange(0.5, len(res.columns), 1), res.columns, rotation=90)
plt.yticks(np.arange(0.5, len(res.index), 1), res.index)
plt.colorbar()
plt.show()

# import seaborn as sns
# import scipy.spatial as sp, scipy.cluster.hierarchy as hc
#
# dist = 1 - res.values * res.values.transpose()
# tmp = pd.DataFrame(dist, index=res.index, columns=res.columns)
# linkage = hc.linkage(sp.distance.squareform(tmp), method='average')
#
# sns.clustermap(tmp, row_linkage=linkage, col_linkage=linkage, cmap="Purples")
# plt.show()

# list(itertools.product(["a", "b", "c", "d", "e"], list(range(1,20))))

# for train_index, test_index in kf.split(X):
#     # print("TRAIN:", train_index, "TEST:", test_index)
#     X_train, X_validation = X[train_index], X[test_index]
#     y_train, y_validation = y[train_index], y[test_index]
#
#     clfLDA = LinearDiscriminantAnalysis()
#     clfLDA.fit(X_train, y_train)
#     predsLDA = clfLDA.predict(X_validation)

    # clfNB = GaussianNB()
    # clfNB.fit(X_train, y_train)
    # predsNB = clfNB.predict(X_validation)
    #
    # clfLR = LogisticRegression(solver="liblinear")
    # clfLR.fit(X, y)
    # predsLR = clfLR.predict(X_validation)
    #
    # print(f"LDA: {np.around((1 - f1_score(y_validation, predsLDA))*100, decimals=1)}\n" +
    #       f" NB: {np.around((1 - f1_score(y_validation, predsNB))*100, decimals=1)}\n" +
    #       f" LR: {np.around((1 - f1_score(y_validation, predsLR))*100, decimals=1)}\n")






