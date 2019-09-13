include: "subworkflows.smk"

from scripts.utils import Trigger

DATASET = "neuropeptides"

split_normalize = Trigger(
    swf=split_normalize_swf,
    output=expand("data/{dataset}_{part}/seqs.fasta", dataset=DATASET, part=["ds1", "ds2"]) + \
          expand("data/{dataset}_{part}/classes.txt", dataset=DATASET, part=["ds1", "ds2"]))

sequence_length_distribution = Trigger(
    swf=sequence_length_distribution_swf,
    output=expand("data/{dataset}/sequence_length_distribution.svg", dataset=DATASET))

DATASET = "neuropeptides_ds3"

multiple_sequence_alignment = Trigger(
    swf=generate_multiple_sequence_alignment_swf,
    output=expand("data/{dataset}/seqs_msa.fasta", dataset=DATASET))

secondary_structure_profile = Trigger(
    swf=secondary_structure_profile_swf,
    output=expand("data/{dataset}/annotated_seqs.fasta", dataset=DATASET) + \
          expand("data/{dataset}/annotated_seqs_msa.fasta", dataset=DATASET) + \
          expand("data/{dataset}/annotated_classes.txt", dataset=DATASET)
)

disorder = Trigger(
    swf=disorder_swf,
    output=expand("data/{dataset}/csv/disorder.csv", dataset=DATASET)
)

rule all:
    input:
         split_normalize.trigger,
         sequence_length_distribution.trigger,
         multiple_sequence_alignment.trigger,
         secondary_structure_profile.trigger,
         disorder.trigger

rule utils_split_normalize:
    input:
        split_normalize.trigger
    output:
        split_normalize.output

rule plots_sequence_length_distribution:
    input:
         rules.utils_split_normalize.output,
         sequence_length_distribution.trigger
    output:
        sequence_length_distribution.output

rule utils_multiple_sequence_alignment:
    input:
        rules.utils_split_normalize.output,
        multiple_sequence_alignment.trigger
    output:
        multiple_sequence_alignment.output

rule utils_secondary_structure_profile:
    input:
        rules.utils_multiple_sequence_alignment.output,
        secondary_structure_profile.trigger
    output:
        secondary_structure_profile.output

rule encodings_disorder:
    input:
        disorder.trigger,
        rules.utils_secondary_structure_profile.output
    output:
        disorder.output

rule generate_dag:
    input:
         rules.plots_sequence_length_distribution.output,
         rules.encodings_disorder.output




