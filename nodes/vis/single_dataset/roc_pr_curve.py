from sklearn.metrics import average_precision_score, roc_curve, roc_auc_score, auc
from sklearn.metrics import precision_recall_curve

import numpy as np
import pandas as pd
import altair as alt

from functools import reduce


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
    df = pd.DataFrame({"x": mean_fpr, "y": mean_tpr, "tprs_lower": tprs_lower, "tprs_upper": tprs_upper})
    df["Encoding"] = encoding
    return df


def is_struc_based(e):
    if "asa" in e:
        return True
    elif "delaunay" in e:
        return True
    elif "disorder" in e:
        return True
    elif "hull" in e:
        return True
    elif "qsar" in e:
        return True
    elif "sse" in e:
        return True
    elif "ta" in e:
        return True
    else:
        return False


def get_data(encodings, dataset):
    test_dfs, prob_dfs, names = [], [], []
    for e in encodings:
        test_dfs += [pd.read_csv(f"data/{dataset}/benchmark/single/y_true_cv_{e}.csv", index_col=0)]
        prob_dfs += [pd.read_csv(f"data/{dataset}/benchmark/single/y_prob_cv_{e}.csv", index_col=0)]
        names += [e]
    df_roc = reduce(lambda a, b: pd.concat([a, b]), map(get_roc_data, test_dfs, prob_dfs, names))
    df_roc["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in df_roc["Encoding"]]
    df_prc = reduce(lambda a, b: pd.concat([a, b]), map(get_pr_data, test_dfs, prob_dfs, names))
    df_prc["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in df_prc["Encoding"]]
    return df_roc, df_prc


def roc_chart(df_f1, dataset):

    indices = df_f1.apply(np.median).sort_values(ascending=False).index.tolist()

    top_six_encodings = indices[:6]
    df_roc_t6, df_prc_t6 = get_data(top_six_encodings, dataset)

    top_three_seq = indices[:3]
    top_three_str = [i for i in indices if is_struc_based(i)][:3]
    df_roc_t33, df_prc_t33 = get_data(top_three_seq + top_three_str, dataset)

    random_guess_line = alt.Chart(
        pd.DataFrame({"x": [0.0, 1.0],
                      "y": [0.0, 1.0]})
    ).mark_line(
        color="lightgrey",
        strokeDash=[3, 1]
    ).encode(
        x="x",
        y="y"
    )

    c1 = alt.Chart(df_roc_t6).mark_line().encode(
        x=alt.X("x:Q", axis=alt.Axis(title=None)),
        y=alt.Y("y:Q", axis=alt.Axis(title="Sensitivity")),
        color=alt.Color("Encoding:N")
    ) + random_guess_line

    c2 = alt.Chart(df_prc_t6).mark_line().encode(
        x=alt.X("x:Q", axis=alt.Axis(title=None)),
        y=alt.Y("y:Q", axis=alt.Axis(title="Precision")),
        color="Encoding:N"
    )

    cA = alt.Chart(df_roc_t33).mark_line().encode(
        x=alt.X("x:Q", axis=alt.Axis(title="1 - Specificity")),
        y=alt.Y("y:Q", axis=alt.Axis(title="Sensitivity")),
        color="Encoding:N",
        strokeDash="type:N"
    ) + random_guess_line

    cB = alt.Chart(df_prc_t33).mark_line().encode(
        x=alt.X("x:Q", axis=alt.Axis(title="Recall")),
        y=alt.Y("y:Q", axis=alt.Axis(title="Precision")),
        color="Encoding:N",
        strokeDash="type:N"
    )

    return (c1 | c2) & (cA | cB)
