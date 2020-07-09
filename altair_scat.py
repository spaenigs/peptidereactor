import pandas as pd
import altair as alt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import cross_val_score, RepeatedStratifiedKFold

dataset = "hiv_protease"

e1 = "ngram_a3_300" # "aaindex_BUNA790102"
df1 = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/sequence_based/y_prob_cv_{e1}.csv", index_col=0)
df1a = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/sequence_based/y_true_cv_{e1}.csv", index_col=0)

e2 = "delaunay_cartesian_product"
df2 = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/structure_based/y_prob_cv_{e2}.csv", index_col=0)
df2a = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/structure_based/y_true_cv_{e2}.csv", index_col=0)

e3 = "aaindex_BUNA790102"
df3 = pd.read_csv(f"data/{dataset}/benchmark/ensemble/seq_vs_str/sequence_based/y_prob_cv_{e3}.csv", index_col=0)

df = pd.DataFrame({"x": df1.values.ravel(), "y": df2.values.ravel(), "z": df3.values.ravel(), "class": df1a.values.ravel()})
df = df.dropna()
X, y = df.iloc[:, :-1].values, df["class"].values


def append_values(y, split_id):
    df = pd.DataFrame({f"split_{split_id}": y}).transpose()
    df.columns = [f"y_{i}" for i in range(df.shape[1])]
    return df


brf = \
    RandomForestClassifier(n_estimators=100, random_state=0)

cv = \
    RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)

df_y_true, df_y_pred, df_y_prob, df_imp = \
    pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

for i, (train_index, test_index) in enumerate(cv.split(X, y)):
    X_train, y_train = X[train_index], y[train_index]
    X_test, y_test = X[test_index], y[test_index]

    df_y_true = pd.concat([df_y_true, append_values(y_test, i)])

    # df_y_true = pd.concat([df_y_true, append_values(y_test, i)])
    brf.fit(X_train, y_train)
    y_pred_class = brf.predict(X_test)
    y_pred_proba = brf.predict_proba(X_test)
    df_y_prob = pd.concat([df_y_prob, append_values(y_pred_proba[:, 1], i)])
    print(f1_score(y_test, y_pred_class))


df_y_true.to_csv("y_true.csv")
df_y_prob.to_csv("y_prob.csv")

c = alt.Chart(df).mark_point(shape="cross", size=1).encode(
    x=alt.X("x", title=f"Predicted probability {e1}"),
    y=alt.Y("y", title=f"Predicted probability {e2}"),
    color="class:N"
)

c.save("scatter.html")



