configfile: "config.yaml"

import sys
sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])

import scripts.utils as utils

workdir: "."

include: "01_preprocessing/a_preprocessing.smk"
include: "01_preprocessing/b_profiles.smk"

include: "02_encoding/a_encode.smk"
include: "02_encoding/b_filter_and_normalize.smk"
include: "02_encoding/c_final_datasets.smk"

# TODO structure based
# TODO check aaindex distance (1 -?)

DATASET = config["dataset"]
PART = config["part"]
NORMALIZE = config["normalize"]
# ENCODINGS = ["disorder", "spinex", "pssm", "psipred"]
ENCODINGS = ["psekraac"]
rule all:
    input:
        # expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" +
        #        "{dataset}_{part}_normalized-{normalized}.txt",
        #        dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=ENCODINGS)
        # expand("00_data/out/{dataset}/plots/{dataset}_length_distribution.svg", dataset=DATASET),
        expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" +
               "geom_median/tsne/normalized-{normalized}/final_datasets.txt",
               dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=ENCODINGS),
        # TODO works only for param_based, psekraac and aaindex encoding:
        expand("00_data/out/{dataset}/plots/{dataset}_{part}_{encoding}_normalized-{normalized}_tsne.svg",
               dataset=DATASET, part=PART, normalized=NORMALIZE, encoding=ENCODINGS),








# subworkflow a_preprocessing_preprocessing:
#     workdir:
#         "."
#     snakefile:
#         "02_preprocessing/a_preprocessing.smk"
#     configfile:
#         "config.yaml"
#
#
# subworkflow a_preprocessing_profiles_and_filter:
#     workdir:
#         "."
#     snakefile:
#         "02_preprocessing/b_profiles.smk"
#     configfile:
#         "config.yaml"
#
#
# subworkflow b_psekraac_encode:
#     workdir:
#         "."
#     snakefile:
#         "03_encoding/psekraac/a_encode.smk"
#     configfile:
#         "config.yaml"
#
#
# subworkflow b_psekraac_filter_and_normalize:
#     workdir:
#         "."
#     snakefile:
#         "03_encoding/psekraac/b_filter_and_normalize.smk"
#     configfile:
#         "config.yaml"
#
#
# subworkflow b_psekraac_final_datasets:
#     workdir:
#         "."
#     snakefile:
#         "03_encoding/psekraac/c_final_datasets.smk"
#     configfile:
#         "config.yaml"
#
#
# def target_files(encoding):
#     files = []
#     if encoding == "psekraac":
#         for type_ in config["psekraac"]["types"]:
#             files += expand("01_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/csv/original/" + \
#                                 "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}.csv",
#                             dataset=config["dataset"], part=config["part"],
#                             name=config["psekraac"][type_]["name"],
#                             subtype=config["psekraac"][type_]["subtypes"],
#                             raactype=config["psekraac"][type_]["raactypes"],
#                             ktuple=config["psekraac"][type_]["ktuples"],
#                             glambda=config["psekraac"][type_]["glambdas"])[:10]
#     print(len(files))
#     return files
#
# rule all:
#     input:
#         a_preprocessing_preprocessing(f"01_data/out/{DATASET}/fasta/{DATASET}_{PART}.fasta"),
#         a_preprocessing_preprocessing(f"01_data/out/{DATASET}/class/{DATASET}_{PART}_classes.txt"),
#         a_preprocessing_preprocessing(f"01_data/out/{DATASET}/plots/{DATASET}_length_distribution.svg"),
#         a_preprocessing_profiles_and_filter(
#             f"01_data/out/{DATASET}/{DATASET}_{PART}/joblib/{DATASET}_{PART}_pssms_filtered.joblib"),
#         a_preprocessing_profiles_and_filter(
#             f"01_data/out/{DATASET}/{DATASET}_{PART}/joblib/{DATASET}_{PART}_pssms_filtered_msa.joblib"),
#         b_psekraac_encode(target_files("psekraac"))
#         b_psekraac_filter_and_normalize(
#             expand("01_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/csv/normalized/{dataset}_{part}_normalized-{normalized}.txt",
#                dataset=DATASET, part=PART, normalized=NORMALIZE)),
#         b_psekraac_final_datasets(
#             expand("01_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}.txt",
#                dataset=DATASET, part=PART, normalized=NORMALIZE)),
#         b_psekraac_final_datasets(expand("01_data/out/{dataset}/plots/{dataset}_{part}_normalized-{normalized}_tsne.svg",
#                                          dataset=DATASET, part=PART, normalized=NORMALIZE))




