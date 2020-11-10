import pandas as pd
import numpy as np
from scipy import stats

DATASETS = {
    "cpp_mixed": -1,
    "amp_gonzales": -1,
    "cpp_sanders": -1,
    "hiv_bevirimat": -1,
    "tce_zhao": -1,
    "amp_fernandes": -1,
    "hiv_nfv": -1,
    "hiv_v3": -1,
    "amp_csamp": -1,
    "acp_iacp": -1,
    "cpp_mlcppue": -1,
    "acp_anticp": -1,
    "atb_antitbp": -1,
    "hiv_lpv": -1,
    "acp_mlacp": -1,
    "hiv_abc": -1,
    "hiv_azt": -1,
    "hiv_d4t": -1,
    "hiv_ddi": -1,
    "hiv_3tc": -1,
    "amp_iamp2l": -1,
    "hem_hemopi": -1,
    "aip_aippred": -1,
    "hiv_apv": -1,
    "hiv_dlv": -1,
    "hiv_efv": -1,
    "hiv_rtv": -1,
    "hiv_nvp": -1,
    "hiv_idv": -1,
    "amp_antibp": -1,
    "cpp_cppredfl": -1,
    "hiv_protease": -1,
    "cpp_kelmcpp": -1,
    "avp_avppred": -1,
    "cpp_cellppdmod": -1,
    "pip_pipel": -1,
    "hiv_sqv": -1,
    "atb_iantitb": -1,
    "avp_amppred": -1,
    "cpp_cellppd": -1,
    "nep_neuropipred": -1,
    "cpp_mlcpp": -1,
    "amp_antibp2": -1,
    "bce_ibce": -1,
    "amp_modlamp": -1,
    "afp_amppred": -1,
    "afp_antifp": -1,
    "isp_il10pred": -1,
    "ace_vaxinpad": 0.95,  # Acc
    "aip_antiinflam": -1
}

for ds in DATASETS.keys():
    df = pd.read_csv(f"data/ds/benchmark/metrics/f1.csv", index_col=0)
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
    print(f"{ds},{DATASETS[ds]},{np.mean(f1s)} (+- {std}),{sig_lvl}")