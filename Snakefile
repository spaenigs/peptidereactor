configfile: "config.yaml"

import sys
sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])

try:
    import encoder.ifeature.param_free.encoder as param_free_encoder
except ModuleNotFoundError as e:
    print()
    print(e)
    print("Loaded environment specific configuration?", file=sys.stderr)
    sys.exit(1)

import scripts.utils as utils

workdir: "."

# include: "01_preprocessing/a_preprocessing.smk"
# include: "01_preprocessing/b_profiles.smk"
#
# include: "02_encoding/a_encode.smk"
# include: "02_encoding/b_normalize.smk"
# include: "02_encoding/c_final_datasets.smk"

include: "03_machine_learning/a_train_test_split.smk"

DATASET = config["dataset"]
PART = config["part"]
NORMALIZE = config["normalize"]

ENCODINGS = sorted(utils.STRUC_ENCODINGS +
                   utils.REST_ENCODINGS +
                   utils.PARAM_BASED_ENCODINGS +
                   utils.PARAM_FREE_ENCODINGS)

ENCODINGS_PLOT = sorted([utils.PSEKRAAC] + utils.PARAM_BASED_ENCODINGS)

# ENCODINGS_TTEST = utils.

# ENCODINGS = utils.STRUC_ENCODINGS + utils.PARAM_FREE_ENCODINGS
# ENCODINGS_PLOT = []

rule all:
    input:
        # expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" +
        #        "{dataset}_{part}_normalized-{normalized}.txt",
        #        dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=ENCODINGS)
        # expand("00_data/out/{dataset}/plots/{dataset}_length_distribution.svg", dataset=DATASET),
        # expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" +
        #        "geom_median/tsne/normalized-{normalized}/final_datasets.txt",
        #        dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=utils.APAAC),
        # # TODO works only for param_based, psekraac and aaindex encoding:
        # expand("00_data/out/{dataset}/plots/{dataset}_{part}_{encoding}_normalized-{normalized}_tsne.svg",
        #        dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=utils.PSEKRAAC),
        # expand("00_data/out/{dataset}/plots/{encoding}/{dataset}_{part}_normalized-yes_ttest_error.pdf",
        #        dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=[utils.APAAC, utils.BINARY, utils.PSEKRAAC]),
        expand("00_data/out/{dataset}/{dataset}_{part}/analyis/t_test/" + \
               "{dataset}_{part}_normalized-{normalized}_ttest_error_all.csv",
               dataset=DATASET, part=PART, normalized=NORMALIZE),
         expand("00_data/out/neuropeptides/plots/" + \
                "{dataset}_{part}_normalized-{normalized}_ttest_error_all_vs_all.pdf",
                dataset=DATASET, part=PART, normalized=NORMALIZE)