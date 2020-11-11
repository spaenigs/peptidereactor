import pandas as pd
import numpy as np
from scipy import stats

DATASETS = {
    "cpp_mixed": -1,  # 0.83,  # Acc, use with caution!
    "amp_gonzales": -1,
    "cpp_sanders": 0.9172,  # Acc
    "hiv_bevirimat": 0.927,  # Acc
    "tce_zhao": 0.833,  # AUC
    "amp_fernandes": 0.967,  # Acc
    "hiv_nfv": 0.96,  # Acc
    "hiv_v3": 0.937,  # AUC
    "amp_csamp": 0.9,  # Acc
    "acp_iacp": 0.9506,  # Acc
    "cpp_mlcppue": 0.725,  # Acc
    "acp_anticp": 0.9144,  # Acc
    "atb_antitbp": 0.7854,  # Acc
    "hiv_lpv": 0.98,  # Acc
    "acp_mlacp": 0.887,  # Acc
    "hiv_abc": 0.95,  # Acc
    "hiv_azt": 0.94,  # Acc
    "hiv_d4t": 0.87,  # Acc
    "hiv_ddi": 0.85,  # Acc
    "hiv_3tc": 0.97,  # Acc
    "amp_iamp2l": 0.9223,  # Acc
    "hem_hemopi": 0.953,  # Acc
    "aip_aippred": 0.814,  # Acc
    "hiv_apv": 0.91,  # Acc
    "hiv_dlv": 0.92,  # Acc
    "hiv_efv": 0.94,  # Acc
    "hiv_rtv": 0.97,  # Acc
    "hiv_nvp": 0.96,  # Acc
    "hiv_idv": 0.97,  # Acc
    "amp_antibp": 0.9211,  # Acc
    "cpp_cppredfl": 0.921,  # Acc
    "hiv_protease": 0.932,  # Acc
    "cpp_kelmcpp": 0.8698,  # Acc
    "avp_avppred": 0.85,  # Acc
    "cpp_cellppdmod": 0.951,  # Acc
    "pip_pipel": 0.717,  # Acc
    "hiv_sqv": 0.9,  # Acc
    "atb_iantitb": 0.862,  # Acc
    "avp_amppred": 0.9009,
    "cpp_cellppd": 0.9740,  # Acc
    "nep_neuropipred": 0.8650,  # Acc, ds1 only
    "cpp_mlcpp": 0.896,  # Acc
    "amp_antibp2": 0.9255,  # Acc
    "bce_ibce": 0.732,  # Acc
    "amp_modlamp": -1,
    "afp_amppred": 0.9335,  # Acc
    "afp_antifp": 0.8878,  # Acc
    "isp_il10pred": 0.8124,  # Acc
    "ace_vaxinpad": 0.95,  # Acc
    "aip_antiinflam": 0.781  # Acc
}

res = []
for ds in DATASETS.keys():
    df = pd.read_csv(f"data/{ds}/benchmark/metrics/f1.csv", index_col=0)
    best_encoding = df.apply(np.median).sort_values(ascending=False).index[0]
    f1s = df.loc[:, best_encoding]
    std = np.std(f1s)
    _, pvalue = stats.ttest_1samp(f1s, DATASETS[ds])
    sig_lvl = \
        "****" if pvalue <= 0.0001 else \
        "***" if pvalue <= 0.001 else \
        "**" if pvalue <= 0.01 else \
        "*" if pvalue <= 0.05 else \
        "ns."
    diff = np.abs(np.mean(f1s) - DATASETS[ds])
    # diff =  DATASETS[ds] - np.mean(f1s)
    res += [[ds, DATASETS[ds], np.mean(f1s), std, sig_lvl, diff]]

df_res = pd.DataFrame(res)
df_res.columns = ["dataset", "performance", "f1", "std", "sig", "diff"]
df_res.sort_values(by="diff", inplace=True)

# top10 = df_res.iloc[:, :].sort_values(by="f1", ascending=False)
top10 = df_res.iloc[:10, :].sort_values(by="f1", ascending=False)
# print(top10)

for n, s in top10.iterrows():
    if s["performance"] == -1.0:
        continue
    ds = s['dataset'].replace('_', '\\_')
    sig_lvl = s["sig"].replace("*", "\\ast")
    perf_tmp = np.round(s["performance"], decimals=3)
    perf = f"\\textbf{{{perf_tmp}}}" if s["performance"] > s["f1"] else perf_tmp
    f1_tmp = np.round(s["f1"], decimals=3)
    f1 = f"\\textbf{{{f1_tmp}}}" if s["f1"] > s["performance"] else f1_tmp
    std_tmp = np.round(s['std'], decimals=2)
    print(f"{ds} & {perf} & {f1} ({std_tmp}) & ${sig_lvl}$ & (66) \\\\")
