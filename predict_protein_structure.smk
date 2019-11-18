import os
from modlamp.core import read_fasta

config["global_workdir"] = os.getcwd() + "/"

DATASET = config["dataset"]

rule all:
    input:
         expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])

rule util_protein_structure_prediction:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         link_in="nodes/utils/protein_structure_prediction/raptorx_download_link.txt"
    output:
         pdbs_out=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1])
    params:
         subworkflow="protein_structure_prediction",
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml"
    resources:
         cores=-1
    script:
         "utils/subworkflow.py"