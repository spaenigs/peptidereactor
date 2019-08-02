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


def compute_all_vs_all_t_test_on_classification_error(file_paths, names):

    res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                       index=names,
                       columns=names)

    dfs = ((pd.read_csv(path, index_col=0), name) for (path, name) in zip(file_paths, names))

    for (df_tmp1, n1), (df_tmp2, n2) in itertools.combinations(dfs, 2):
        errors_e1, N_training, N_testing = get_classification_error(df_tmp1)
        errors_e2, _, _ = get_classification_error(df_tmp2)
        m_diff = np.abs(np.mean(errors_e1 - errors_e2))
        sd = np.std(errors_e1 - errors_e2, ddof=1)
        standard_error_mean_corr = sd * np.sqrt((1 / len(errors_e1)) + (N_testing / N_training))
        df = len(errors_e1) - 1
        res.loc[n1, n2] = 2 * stats.t.cdf(-m_diff / standard_error_mean_corr, df)

    # https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
    fig, ax = plt.subplots()
    pc = ax.imshow(res.values * res.values.transpose(), vmin=0, vmax=1, cmap='Greys')
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

    plt.title("Amphiphilic Pseudo-Amino Acid Composition (APAAC):\n" +
              "paired t-test with corrected variance on classification error")

    plt.show()


PARAM_BASED_ENCODINGS = ["apaac", "paac", "cksaagp", "cksaap", "ctriad",
                         "ksctriad", "geary", "moran", "nmbroto", "qsorder",
                         "socnumber", "eaac", "egaac"]

PARAM_FREE_ENCODINGS = ["binary", "aac", "gaac", "ctdt", "ctdc", "ctdd", "tpc",
                        "gtpc", "gdpc", "dpc", "gdpc", "dde", "blosum62", "zscale"]

REST_ENCODINGS = ["aaindex", "psekraac"]

STRUC_ENCODINGS = ["disorder", "spinex", "psipred", "pssm"]


# file_paths = sorted(glob.glob("00_data/out/neuropeptides/neuropeptides_ds1/encodings/apaac/csv/normalized/*.csv"),
#                         key=lambda x: int(re.match(".*?(\d{1,2}).csv", x).group(1)))
#
# names = list(map(lambda p: re.match(".*?apaacencoder_(.*).csv", p).group(1), file_paths))

# file_paths_aaindex = sorted(glob.glob("00_data/out/neuropeptides/neuropeptides_ds1/encodings/aaindex/csv/normalized/*.csv"))
# names_aaindex = list(map(lambda p: re.match(".*?aaindexencoder_aaindex-(.*?\d+).csv", p).group(1), file_paths_aaindex))
# print(names_aaindex)
# compute_all_vs_all_t_test_on_classification_error(file_paths_aaindex, names_aaindex)

file_paths_psekraac = sorted(glob.glob("00_data/out/neuropeptides/neuropeptides_ds1/encodings/psekraac/csv/final/geom_median/tsne/normalized-yes/*.csv"))
names_psekraac = list(map(lambda p: re.match(".*?psekraac(.*?)_subtype.*", p).group(1), file_paths_psekraac))

