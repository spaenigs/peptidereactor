import os
from modlamp.core import read_fasta

config["global_workdir"] = os.getcwd() + "/"

DATASET = config["dataset"]

rule all:
    input:
         f"data/{DATASET}/csv/disorder.csv"

rule util_multiple_sequence_alignment:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt"
    output:
         fasta_out=f"data/{DATASET}/seqs_msa.fasta"
    params:
         subworkflow="multiple_sequence_alignment",
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    script:
         "utils/subworkflow.py"

rule util_secondary_structure:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         fasta_msa_in=f"data/{DATASET}/seqs_msa.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         # refseq_db="/media/spaenigs/4B1DB7375F3291A11/uniprot_db/uniref90/uniref90.db"
         refseq_db="db_tmp/uniref90.db"
    output:
         fasta_anno_out=f"data/{DATASET}/annotated_seqs.fasta",
         fasta_anno_msa_out=f"data/{DATASET}/annotated_seqs_msa.fasta",
         classes_anno=f"data/{DATASET}/annotated_classes.txt",
         profiles_out=expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
                         ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"])
    params:
         subworkflow="secondary_structure_profile",
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    resources:
         cores=-1
    script:
         "utils/subworkflow.py"

rule encoding_disorder:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/disorder.csv"
    params:
         subworkflow="disorder",
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml",
    script:
         "utils/subworkflow.py"