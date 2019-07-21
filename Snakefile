import os

from Bio import SeqIO
from snakemake.io import expand

configfile: "config.yaml"

include: "02_preprocessing/a_preprocessing.smk"
include: "02_preprocessing/b_profiles.sf"

include: "03_encoding/psekraac/a_encode.sf"
include: "03_encoding/psekraac/b_filter_and_normalize.sf"
include: "03_encoding/psekraac/c_final_datasets.sf"

DATASET = config["dataset"]
PART = config["part"]
NORMALIZE = config["normalize"]


rule all:
    input:
        expand("01_data/out/{dataset}/plots/{dataset}_length_distribution.svg", dataset=DATASET),
        expand("01_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}.txt",
               dataset=DATASET, part=PART, normalized=NORMALIZE),
        expand("01_data/out/{dataset}/plots/{dataset}_{part}_normalized-{normalized}_tsne.svg",
               dataset=DATASET, part=PART, normalized=NORMALIZE)


# subworkflow a_preprocessing_preprocessing:
#     workdir:
#         "."
#     snakefile:
#         "02_preprocessing/01_preprocessing.sf"
#     configfile:
#         "config.yaml"
#
#
# subworkflow a_preprocessing_profiles_and_filter:
#     workdir:
#         "."
#     snakefile:
#         "02_preprocessing/02_profiles.sf"
#     configfile:
#         "config.yaml"
#
#
# subworkflow b_psekraac_encode:
#     workdir:
#         "."
#     snakefile:
#         "03_encoding/psekraac/01_psekraac.sf"
#     configfile:
#         "config.yaml"
#
#
# subworkflow b_psekraac_filter_and_normalize:
#     workdir:
#         "."
#     snakefile:
#         "03_encoding/psekraac/02_filter_and_normalize.sf"
#     configfile:
#         "config.yaml"
#
#
# subworkflow b_psekraac_final_datasets:
#     workdir:
#         "."
#     snakefile:
#         "03_encoding/psekraac/03_final_datasets.sf"
#     configfile:
#         "config.yaml"


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
#                             glambda=config["psekraac"][type_]["glambdas"])
#     return files



# rule all:
#     input:
#         a_preprocessing_preprocessing(f"01_data/out/{DATASET}/fasta/{DATASET}_{PART}.fasta"),
#         a_preprocessing_preprocessing(f"01_data/out/{DATASET}/class/{DATASET}_{PART}_classes.txt"),
#         a_preprocessing_preprocessing(f"01_data/out/{DATASET}/plots/{DATASET}_length_distribution.svg"),
#
#         a_preprocessing_profiles_and_filter(
#             f"01_data/out/{DATASET}/{DATASET}_{PART}/joblib/{DATASET}_{PART}_pssms_filtered.joblib"),
#         a_preprocessing_profiles_and_filter(
#             f"01_data/out/{DATASET}/{DATASET}_{PART}/joblib/{DATASET}_{PART}_pssms_filtered_msa.joblib"),
#
#         b_psekraac_encode(target_files("psekraac")),
#         b_psekraac_filter_and_normalize(
#             expand("01_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/csv/normalized/{dataset}_{part}_normalized-{normalized}.txt",
#                dataset=DATASET, part=PART, normalized=NORMALIZE)),
#         b_psekraac_final_datasets(
#             expand("01_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}.txt",
#                dataset=DATASET, part=PART, normalized=NORMALIZE)),
#         b_psekraac_final_datasets(expand("01_data/out/{dataset}/plots/{dataset}_{part}_normalized-{normalized}_tsne.svg",
#                                          dataset=DATASET, part=PART, normalized=NORMALIZE))


