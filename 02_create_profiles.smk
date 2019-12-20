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
         expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
                seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1],
                ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"]),
         f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         f"data/{DATASET}/annotated_pdbs_classes.txt",
         expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])

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
             shell(f"""snakemake -s {{params.snakefile}} {{output.fasta_out}} \
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
         profiles_out=expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
                         ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"])
    params:
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.fasta_anno_out}} {{output.fasta_anno_msa_out}} {{output.classes_anno}} {{output.profiles_out}} \
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
         pdbs_out=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1])
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.fasta_out}} {{output.classes_out}} {{output.pdbs_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")