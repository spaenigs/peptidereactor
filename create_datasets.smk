import os

config["global_workdir"] = os.getcwd() + "/"

rule all:
    input:
         "data/neuropeptides/sequence_length_distribution.svg"

rule plot_sequence_length_distribution:
        input:
             fasta_in="data/neuropeptides/seqs.fasta",
             classes_in="data/neuropeptides/classes.txt",
             fasta_out_1="data/neuropeptides_ds1/seqs.fasta",
             classes_out_1="data/neuropeptides_ds1/classes.txt",
             fasta_out_2="data/neuropeptides_ds2/seqs.fasta",
             classes_out_2="data/neuropeptides_ds2/classes.txt"
        output:
             svg_out="data/neuropeptides/sequence_length_distribution.svg"
        params:
             subworkflow="sequence_length_distribution",
             snakefile="nodes/plots/sequence_length_distribution/Snakefile",
             configfile="nodes/plots/sequence_length_distribution/config.yaml"
        script:
             "utils/subworkflow.py"

rule util_split_normalize:
    input:
         fasta_in="data/neuropeptides/seqs.fasta",
         classes_in="data/neuropeptides/classes.txt"
    output:
         fasta_out_1="data/neuropeptides_ds1/seqs.fasta",
         classes_out_1="data/neuropeptides_ds1/classes.txt",
         fasta_out_2="data/neuropeptides_ds2/seqs.fasta",
         classes_out_2="data/neuropeptides_ds2/classes.txt"
    params:
         subworkflow="split_normalize",
         snakefile="nodes/utils/split_normalize/Snakefile",
         configfile="nodes/utils/split_normalize/config.yaml"
    script:
         "utils/subworkflow.py"