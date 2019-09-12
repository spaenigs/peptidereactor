class Trigger:

    def trigger(self, wildcards):
        return self.swf(self.files)

    def __init__(self, swf, files):
        self.files = files
        self.swf = swf

subworkflow split_normalize_swf:
    workdir:
        "."
    snakefile:
        "nodes/utils/split_normalize/Snakefile"
    # configfile:
    #     "nodes/utils/split_normalize/config.yaml"

subworkflow sequence_length_distribution_swf:
    workdir:
        "."
    snakefile:
        "nodes/plots/sequence_length_distribution/Snakefile"
    # configfile:
    #     "nodes/plots/sequence_length_distribution/config.yaml"

DATASET = "neuropeptides"

sequence_length_distribution = Trigger(
    sequence_length_distribution_swf,
    expand("data/{dataset}/sequence_length_distribution.svg", dataset=DATASET))

split_normalize = Trigger(
    split_normalize_swf,
    expand("data/{dataset}_{part}/seqs.fasta", dataset=DATASET, part=["ds1", "ds2"]) + \
    expand("data/{dataset}_{part}/classes.txt", dataset=DATASET, part=["ds1", "ds2"]))

rule all:
    input:
         split_normalize.trigger,
         sequence_length_distribution.trigger

rule utils_split_normalize:
    input:
        split_normalize.trigger
    output:
        split_normalize.files

rule plots_sequence_length_distribution:
    input:
         rules.utils_split_normalize.output,
         sequence_length_distribution.trigger
    output:
        sequence_length_distribution.files

rule generate_dag:
    input:
         rules.plots_sequence_length_distribution.output




