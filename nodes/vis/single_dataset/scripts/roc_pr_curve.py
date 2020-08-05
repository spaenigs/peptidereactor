from sklearn.metrics import *

import numpy as np
import pandas as pd


def get_pr_data(df_test, df_prob, encoding):
    precs = []
    aps = []
    mean_recall = np.linspace(0, 1, 50)
    for i in range(df_test.shape[0]):
        y_true = df_test.iloc[i, :].dropna().values
        y_pred = df_prob.iloc[i, :].dropna().values
        prec, recall, _ = precision_recall_curve(y_true, y_pred)
        interp_prec = np.interp(mean_recall, list(reversed(recall)), list(reversed(prec)))
        precs.append(interp_prec)
        aps.append(average_precision_score(y_true, y_pred))
    mean_precs = np.mean(precs, axis=0)
    df = pd.DataFrame({"x": mean_recall, "y": mean_precs})
    df["Encoding"] = encoding
    df["mean_ap"] = np.round(np.mean(aps), 2)
    df["legend_label"] = df.apply(lambda row: f"{row['Encoding']} (AP: {row['mean_ap']})", axis=1)
    return df


def get_roc_data(df_test, df_prob, encoding):
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    for i in range(df_test.shape[0]):
        y_true = df_test.iloc[i, :].dropna().values
        y_pred = df_prob.iloc[i, :].dropna().values
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(roc_auc_score(y_true, y_pred))
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    df = pd.DataFrame({
        "x": mean_fpr,
        "y": mean_tpr,
        "tprs_lower": tprs_lower,
        "tprs_upper": tprs_upper
    })
    df["Encoding"] = encoding
    df["mean_auc"] = np.round(mean_auc, 2)
    df["legend_label"] = df.apply(lambda row: f"{row['Encoding']} (AUC: {row['mean_auc']})", axis=1)
    return df
