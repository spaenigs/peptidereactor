import os

config["dataset"] = "neuropeptides_ds3"
config["global_workdir"] = os.getcwd() + "/"

rule all:
    input:
        "data/neuropeptides_ds3/csv/disorder.csv",
        "data/neuropeptides_ds3/csv/disorderb.csv",
        "data/neuropeptides_ds3/csv/disorderc.csv",
        # "data/neuropeptides_ds3/csv/aac.csv"

rule encoding_disorder:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorder.csv"
    params:
         subworkflow="disorder",
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_disorderb:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs_msa.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorderb.csv"
    params:
         subworkflow="disorderb",
         snakefile="nodes/encodings/disorderb/Snakefile",
         configfile="nodes/encodings/disorderb/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_disorderc:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorderc.csv"
    params:
         subworkflow="disorderc",
         snakefile="nodes/encodings/disorderc/Snakefile",
         configfile="nodes/encodings/disorderc/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_aac:
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
