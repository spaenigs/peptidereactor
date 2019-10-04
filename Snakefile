# include: "subworkflows.smk"

from utils.utils import Trigger

DATASET = "neuropeptides"

# split_normalize = Trigger(
#     swf=split_normalize_swf,
#     output=expand("data/{dataset}_{part}/seqs.fasta", dataset=DATASET, part=["ds1", "ds2"]) + \
#           expand("data/{dataset}_{part}/classes.txt", dataset=DATASET, part=["ds1", "ds2"]))
#
# sequence_length_distribution = Trigger(
#     swf=sequence_length_distribution_swf,
#     output=expand("data/{dataset}/sequence_length_distribution.svg", dataset=DATASET))
#
DATASET = "neuropeptides_ds3"
#
# multiple_sequence_alignment = Trigger(
#     swf=generate_multiple_sequence_alignment_swf,
#     output=expand("data/{dataset}/seqs_msa.fasta", dataset=DATASET))
#
# secondary_structure_profile = Trigger(
#     swf=secondary_structure_profile_swf,
#     output=expand("data/{dataset}/annotated_seqs.fasta", dataset=DATASET) + \
#           expand("data/{dataset}/annotated_seqs_msa.fasta", dataset=DATASET) + \
#           expand("data/{dataset}/annotated_classes.txt", dataset=DATASET)
# )
#
# disorder = Trigger(
#     swf=disorder_swf,
#     output=expand("data/{dataset}/csv/disorder.csv", dataset=DATASET)
# )
#
# disorderb = Trigger(
#     swf=disorderb_swf,
#     output=expand("data/{dataset}/csv/disorderb.csv", dataset=DATASET)
# )
#
# disorderc = Trigger(
#     swf=disorderc_swf,
#     output=expand("data/{dataset}/csv/disorderc.csv", dataset=DATASET)
# )
#
# aac = Trigger(
#     swf=aac_swf,
#     output=expand("data/{dataset}/csv/aac.csv", dataset=DATASET)
# )

## validation of input and outputfiles on-the-fly

# from snakemake.workflow import Subworkflow
# swf = Subworkflow("all", "disorderc_swf","nodes/encodings/disorderc/Snakefile",
#                   ".", "nodes/encodings/aac/config.yaml")
# print(swf.target())
# import inspect
# print(inspect(disorderc_swf))
# print(inspect.getmembers(aac_swf)) # , predicate=inspect.isclass
# print(aac_swf.__self__.configfile)
# print(aac_swf.__self__.snakefile)
# print(inspect.getmembers(aac_swf.__self__))

# with open(aac_swf.__self__.snakefile) as f:
#     print(f.readlines())

# import yaml

# check how to replace input and output files on runtime...
# with open(aac_swf.__self__.configfile) as f:
#     print(yaml.safe_load(f))

subworkflow aac_swf:
    workdir:
           "."
    snakefile:
           "nodes/encodings/aac/Snakefile"
    configfile:
           "nodes/encodings/aac/config.yaml"

rule all:
    input:
        aac_swf(expand("data/{dataset}/csv/aac.csv", dataset=DATASET))
        # aac.trigger()
         # split_normalize.trigger,
         # sequence_length_distribution.trigger,
         # multiple_sequence_alignment.trigger,
         # secondary_structure_profile.trigger,
         # disorder.trigger,
         # disorderb.trigger,
         # disorderc.trigger
#
# rule utils_split_normalize:
#     input:
#         split_normalize.trigger
#     output:
#         split_normalize.output
#
# rule plots_sequence_length_distribution:
#     input:
#          rules.utils_split_normalize.output,
#          sequence_length_distribution.trigger
#     output:
#         sequence_length_distribution.output
#
# rule utils_multiple_sequence_alignment:
#     input:
#         rules.utils_split_normalize.output,
#         multiple_sequence_alignment.trigger
#     output:
#         multiple_sequence_alignment.output
#
# rule utils_secondary_structure_profile:
#     input:
#         rules.utils_multiple_sequence_alignment.output,
#         secondary_structure_profile.trigger
#     output:
#         secondary_structure_profile.output
#
# rule encodings_disorder:
#     input:
#         disorder.trigger,
#         rules.utils_secondary_structure_profile.output
#     output:
#         disorder.output
#
# rule encodings_disorderb:
#     input:
#         disorderb.trigger,
#         rules.utils_secondary_structure_profile.output
#     output:
#         disorderb.output
#
# rule encodings_disorderc:
#     input:
#         disorderc.trigger,
#         rules.utils_secondary_structure_profile.output
#     output:
#         disorderc.output
#
# rule encoding_benchmark:
#     input:
#          rules.plots_sequence_length_distribution.output,
#          rules.encodings_disorder.output,
#          rules.encodings_disorderb.output,
#          rules.encodings_disorderc.output
#
#
#
#
