import os

config["global_workdir"] = os.getcwd() + "/"

rule all:
    input:
        "data/neuropeptides_ds3/csv/disorder.csv",
        "data/neuropeptides_ds3/csv/aac.csv"

rule a:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorder.csv"
    params:
         subworkflow="disorder",
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml"
    script:
          "utils/subworkflow.py"

rule b:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/aac.csv"
    params:
         subworkflow="aac",
         snakefile="nodes/encodings/aac/Snakefile",
         configfile="nodes/encodings/aac/config.yaml"
    script:
        "utils/subworkflow.py"
