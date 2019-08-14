import scripts.utils as utils
import pandas as pd


rule t_test_on_classification_error:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest.csv"
    script:
        "scripts/i_t_test_on_classification_error.py"


rule plot_t_test_on_classification_error:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error.csv"
    output:
        "00_data/out/neuropeptides/plots/{encoding}/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest.pdf"
    script:
          "scripts/ii_plot_t_test.py"


rule combine_cv:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/{computation}/" + \
                                 "{dataset}_{part}_normalized-{normalized}_cross_validation.csv",
                                 dataset=wildcards.dataset,
                                 part=wildcards.part,
                                 # TODO
                                 encoding=utils.PARAM_FREE_ENCODINGS + \
                                          utils.REST_ENCODINGS + utils.STRUC_ENCODINGS + \
                                          ["apaac", "paac",
                                           "geary", "moran", "nmbroto", "qsorder",
                                           "socnumber"],
                                 computation=wildcards.computation,
                                 normalized=wildcards.normalized)
    output:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    run:
        res = pd.DataFrame()
        for path in list(input):
            res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0, sort=True)
        res.to_csv(str(output))


rule all_vs_all_ttest:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_cross_validation_all.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/{computation}/" + \
         "{dataset}_{part}_normalized-{normalized}_ttest_all.csv"
    threads:
        8
    script:
        "scripts/iii_all_vs_all_ttest.py"


rule generate_network_data:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_all.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/{computation}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_all.json"
    script:
         "scripts/iv_generate_network_data.py"
