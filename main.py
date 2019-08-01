import pandas as pd
import numpy as np

from scipy import stats

import itertools
import re
import glob

from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.metrics import f1_score
import matplotlib.pyplot as plt

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


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


def compute_all_vs_all_t_test_on_classification_error(file_paths):

    def get_dataframes(paths):
        for path in paths:
            yield pd.read_csv(path, index_col=0), path.replace("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/neuropeptides_ds1_apaacencoder_", "").replace(".csv", "")

    idx = list(map(lambda p: p.replace("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/neuropeptides_ds1_apaacencoder_", "").replace(".csv", ""), file_paths))

    res = pd.DataFrame(np.zeros((30, 30)) + 1.0,
                       index=idx,
                       columns=idx)

    for (df_tmp1, n1), (df_tmp2, n2) in itertools.combinations(get_dataframes(file_paths), 2):
        errors_e1, N_training, N_testing = get_classification_error(df_tmp1)
        errors_e2, _, _ = get_classification_error(df_tmp2)
        m_diff = np.abs(np.mean(errors_e1 - errors_e2))
        sd = np.std(errors_e1 - errors_e2, ddof=1)
        standard_error_mean_corr = sd * np.sqrt((1 / len(errors_e1)) + (N_testing / N_training))
        df = len(errors_e1) - 1
        res.loc[n1, n2] = 2 * stats.t.cdf(-m_diff / standard_error_mean_corr, df)

    fig, ax = plt.subplots()
    pc = ax.imshow(res.values * res.values.transpose(), vmin=0, vmax=1, cmap='Greys')
    ax.set_xticks(np.arange(len(res.columns)))
    ax.set_yticks(np.arange(len(res.index)))
    ax.set_xticklabels(res.columns)
    ax.set_yticklabels(res.index)

    # ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    cbar = ax.figure.colorbar(pc)
    cbar.ax.set_ylabel("p-value, (*): <= 0.05, (**): <= 0.01", rotation=-90, va="bottom")

    for i in range(len(res.index)):
        for j in range(len(res.columns)):
            if res.loc[res.index[i], res.columns[j]] <= 0.01:
                # ax.annotate(str(np.round(res.loc[res.index[i], res.columns[j]], decimals=2)),
                #             xy=(j, i), xycoords='data',
                #             xytext=(j + 3, 32), textcoords='data',
                #             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
                #             )
                ax.text(j, i, "(**)", ha="center", va="center", color="black")
            elif res.loc[res.index[i], res.columns[j]] <= 0.05:
                ax.text(j, i, "(*)", ha="center", va="center", color="black")


    fig.tight_layout()
    fig.set_size_inches(13, 6)

    plt.title("Amphiphilic Pseudo-Amino Acid Composition (APAAC)\nPaired t-test with corrected variance on classification error")

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


file_paths = sorted(glob.glob("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/*.csv"),
                        key=lambda x: int(re.match(".*?(\d{1,2}).csv", x).group(1)))

compute_all_vs_all_t_test_on_classification_error(file_paths)


