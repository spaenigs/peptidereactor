from sklearn.metrics import davies_bouldin_score
from glob import glob

import pandas as pd
import numpy as np


def _ravel_and_annotate(df1, df2, df1_class, cat, div, e1, e2, f1_e1, f1_e2):
    df = pd.DataFrame({
        "x": df1.values.ravel(),
        "y": df2.values.ravel(),
        "class": df1_class.values.ravel()})
    df.dropna(inplace=True)
    dbs = davies_bouldin_score(df.iloc[:, :2].values, df["class"])
    df["diversity"] = \
        f"{cat} ({np.round(div, 2)}), " \
        f"e1: {e1} ({np.round(f1_e1, 1)}), " \
        f"e2: {e2} ({np.round(f1_e2, 1)}) - " \
        f"DBS: {np.round(dbs, 2)}"
    df["diversity_e1"], df["diversity_e2"] = f"{e1} - F1={np.round(f1_e1, 1)}", f"{e2} - F1={np.round(f1_e2, 1)}"
    df["diversity_score"], df["diversity_dbs"] = np.round(div, 2), np.round(dbs, 2)
    return df


def _get_data(base_dir, e1, e2):
    df1 = pd.read_csv(glob(base_dir + f"*/y_prob_cv_{e1}.csv")[0], index_col=0)
    df1_class = pd.read_csv(glob(base_dir + f"*/y_true_cv_{e1}.csv")[0], index_col=0)
    df2 = pd.read_csv(glob(base_dir + f"*/y_prob_cv_{e2}.csv")[0], index_col=0)
    return df1, df2, df1_class


def get_crap_combination(base_dir, df_div, medians):
    i_crap, c_crap, div_crap = \
             sorted([(s, i, df_div.loc[s, i])
                     for i, s in df_div.idxmin().items()],
                    key=lambda x: x[2],
                    reverse=True)[0]
    df1_crap, df2_crap, df1_crap_class = _get_data(base_dir, i_crap, c_crap)
    return _ravel_and_annotate(df1_crap, df2_crap, df1_crap_class,
                               cat=f"(a) crap", div=div_crap, e1=i_crap, e2=c_crap,
                               f1_e1=medians[i_crap], f1_e2=medians[c_crap])


def get_low_combination(base_dir, df_div_sub, medians):
    i_low, c_low, div_low = \
        sorted([(s, i, df_div_sub.loc[s, i])
                for i, s in df_div_sub.idxmin().items()],
               key=lambda x: x[2])[0]
    df1_low, df2_low, df1_low_class = _get_data(base_dir, i_low, c_low)
    return _ravel_and_annotate(df1_low, df2_low, df1_low_class,
                               cat=f"(b) low", div=div_low, e1=i_low, e2=c_low,
                               f1_e1=medians[i_low], f1_e2=medians[c_low])


def get_mid_combination(base_dir, df_div_sub, medians):
    s = pd.Series([df_div_sub.loc[i, j]
                   for i, j in [(i, j)
                                for i in df_div_sub.index
                                for j in df_div_sub.columns]])
    min_limit, max_limit = s.describe(percentiles=[0.3, 0.7])[["30%", "70%"]]
    i_mid, c_mid, div_mid = sorted([(i, j, df_div_sub.loc[i, j])
                                    for i, j in [(i, j) for i in df_div_sub.index for j in df_div_sub.columns]
                                    if min_limit < df_div_sub.loc[i, j] < max_limit], key=lambda x: x[2])[0]
    df1_mid, df2_mid, df1_mid_class = _get_data(base_dir, i_mid, c_mid)
    return _ravel_and_annotate(df1_mid, df2_mid, df1_mid_class,
                               cat=f"(c) mid", div=div_mid, e1=i_mid, e2=c_mid,
                               f1_e1=medians[i_mid], f1_e2=medians[c_mid])


def get_high_combination(base_dir, df_div_sub, medians):
    i_high, c_high, div_high = \
        sorted([(s, i, df_div_sub.loc[s, i])
                for i, s in df_div_sub.idxmax().items()],
               key=lambda x: x[2])[-1]
    df1_high, df2_high, df1_high_class = _get_data(base_dir, i_high, c_high)
    return _ravel_and_annotate(df1_high, df2_high, df1_high_class,
                               cat=f"(d) high", div=div_high, e1=i_high, e2=c_high,
                               f1_e1=medians[i_high], f1_e2=medians[c_high])


def get_highest_f1_combination(base_dir, i_highest_f1, c_highest_f1, df_div_sub, medians):
    div_highest_f1 = df_div_sub.loc[i_highest_f1, c_highest_f1]
    df1_highest_f1, df2_highest_f1, df1_highest_f1_class = \
        _get_data(base_dir, i_highest_f1, c_highest_f1)
    return _ravel_and_annotate(df1_highest_f1, df2_highest_f1, df1_highest_f1_class,
                               cat=f"(e) highest f1", div=div_highest_f1, e1=i_highest_f1, e2=c_highest_f1,
                               f1_e1=medians[i_highest_f1], f1_e2=medians[c_highest_f1])
