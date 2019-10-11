import os

config["global_workdir"] = os.getcwd() + "/"

rule all:
    input:
         expand("data/{dataset}/annotated_seqs.fasta", dataset=config["dataset"]),
         expand("data/{dataset}/annotated_seqs_msa.fasta", dataset=config["dataset"]),
         expand("data/{dataset}/annotated_classes.txt", dataset=config["dataset"]),
         expand("data/{dataset}/profile/profiles.txt", dataset=config["dataset"])

rule util_multiple_sequence_alignment:
    input:
         fasta_in="data/{dataset}/seqs.fasta",
         classes_in="data/{dataset}/classes.txt"
    output:
         fasta_out="data/{dataset}/seqs_msa.fasta"
    params:
         subworkflow="multiple_sequence_alignment",
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    script:
         "utils/subworkflow.py"

rule util_secondary_structure:
    input:
         fasta_in="data/{dataset}/seqs.fasta",
         fasta_msa_in="data/{dataset}/seqs_msa.fasta",
         classes_in="data/{dataset}/classes.txt",
         refseq_db="/media/spaenigs/4B1DB7375F3291A1/uniprot_db/uniref90/uniref90.db"
    output:
         fasta_anno_out="data/{dataset}/annotated_seqs.fasta",
         fasta_anno_msa_out="data/{dataset}/annotated_seqs_msa.fasta",
         classes_anno="data/{dataset}/annotated_classes.txt",
         profile="data/{dataset}/profile/profiles.txt"
    params:
         subworkflow="secondary_structure_profile",
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    resources:
         cores=-1
    script:
         "utils/subworkflow.py"