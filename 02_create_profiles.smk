from utils.snakemake_config import WorkflowExecuter
from modlamp.core import read_fasta

DATASET = config["dataset"]
CORES = int(config["cores"])

rule all:
    input:
         f"data/{DATASET}/seqs_msa.fasta",
         f"data/{DATASET}/annotated_seqs.fasta",
         f"data/{DATASET}/annotated_seqs_msa.fasta",
         f"data/{DATASET}/annotated_classes.txt",
         f"data/{DATASET}/profile/",
         f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         f"data/{DATASET}/annotated_pdbs_classes.txt",
         f"data/{DATASET}/pdb/"

rule util_multiple_sequence_alignment:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt"
    output:
         fasta_out=f"data/{DATASET}/seqs_msa.fasta"
    params:
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} \
                                 --cores {CORES} \
                                 --directory $PWD \
                                 --configfile {{params.configfile}}""")

rule util_secondary_structure_profile:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         fasta_msa_in=f"data/{DATASET}/seqs_msa.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         uniprot90_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                    "download_links/uniprot90_download_link.txt",
         psipred_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                  "download_links/psipred_download_link.txt",
         spineXpublic_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                       "download_links/spineXpublic_download_link.txt",
         VSL2_download_link_in="nodes/utils/secondary_structure_profile/" + \
                               "download_links/VSL2_download_link.txt"
    output:
         fasta_anno_out=f"data/{DATASET}/annotated_seqs.fasta",
         fasta_anno_msa_out=f"data/{DATASET}/annotated_seqs_msa.fasta",
         classes_anno=f"data/{DATASET}/annotated_classes.txt",
         profiles_out=directory(f"data/{DATASET}/profile/")
    params:
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} \
                                 --cores {CORES} \
                                 --directory $PWD \ 
                                 --configfile {{params.configfile}}""")

rule util_protein_structure_prediction:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         download_link_in="nodes/utils/protein_structure_prediction/" + \
                          "download_links/raptorx_download_link.txt",
         license_key_in="nodes/utils/protein_structure_prediction/" + \
                        "download_links/modeller_license_key.txt"
    output:
         fasta_out=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_out=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdbs_out=directory(f"data/{DATASET}/pdb/")
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} \
                                 --cores {CORES} \
                                 --directory $PWD \
                                 --configfile {{params.configfile}}""")