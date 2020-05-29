import pandas as pd
from sklearn.metrics import matthews_corrcoef, f1_score

encodings = ["data/hiv_protease/benchmark/single/y_true_cv_waac_aaindex_QIAN880101.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_t9_st-g-gap_rt-9_ktu-1_la-2.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_t1_st-g-gap_rt-19_ktu-1_la-3.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_socnumber_nlag_5.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_qsorder_nlag_4.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_paac_lambda_3.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_ngram_e3_200.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_ngram_a2_100.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_gdpc.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_flgc_aaindex_FINA910104.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_fft_aaindex_FASG760103.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_egaac_window_29.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_dist_freq_dn_5_dc_20.csv",
             "data/hiv_protease/benchmark/single/y_true_cv_ssec.csv"]

df_res = pd.DataFrame()

for p in encodings:
    df_true = pd.read_csv(p, index_col=0)
    df_pred = pd.read_csv(p.replace("_true_", "_pred_"), index_col=0)
    f1s = []
    for i in range(50):
        y_true = df_true.iloc[i, :].dropna()
        y_pred = df_pred.iloc[i, :].dropna()
        f1s += [f1_score(y_true, y_pred)]
    n = p.replace("data/hiv_protease/benchmark/single/", "").replace(".csv", "").replace("y_true_cv_", "")
    df_res = pd.concat([df_res, pd.DataFrame({n: f1s}).transpose()])

df_res.transpose().to_csv("res.csv")

