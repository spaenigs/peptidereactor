import os
from modlamp.core import read_fasta

config["global_workdir"] = os.getcwd() + "/"

DATASET = config["dataset"]

rule all:
    input:
         f"data/{DATASET}/csv/qsar.csv"
         # expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
         #        seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
         #        ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"])
         # expand(f"data/{DATASET}/pdb/{{seq_name}}_minimized.pdb",
         #        seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[:10]
         # expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
         #        seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[-1]

# rule util_secondary_structure_profile:
#     input:
#          fasta_in=f"data/{DATASET}/seqs.fasta",
#          fasta_msa_in=f"data/{DATASET}/seqs_msa.fasta",
#          classes_in=f"data/{DATASET}/classes.txt",
#          uniprot90_download_link_in="nodes/utils/secondary_structure_profile/" + \
#                                     "download_links/uniprot90_download_link.txt",
#          psipred_download_link_in="nodes/utils/secondary_structure_profile/" + \
#                                   "download_links/psipred_download_link.txt",
#          spineXpublic_download_link_in="nodes/utils/secondary_structure_profile/" + \
#                                        "download_links/spineXpublic_download_link.txt",
#          VSL2_download_link_in="nodes/utils/secondary_structure_profile/" + \
#                                "download_links/VSL2_download_link.txt"
#     output:
#          fasta_anno_out=f"data/{DATASET}/annotated_seqs.fasta",
#          fasta_anno_msa_out=f"data/{DATASET}/annotated_seqs_msa.fasta",
#          classes_anno=f"data/{DATASET}/annotated_classes.txt",
#          profiles_out=expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
#                          seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
#                          ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"])
#     params:
#          subworkflow="secondary_structure_profile",
#          snakefile="nodes/utils/secondary_structure_profile/Snakefile",
#          configfile="nodes/utils/secondary_structure_profile/config.yaml",
#          dryrun=config["dryrun"]
#     resources:
#          cores=1
#     script:
#          "utils/subworkflow.py"

rule util_protein_structure_prediction:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         download_link_in="nodes/utils/protein_structure_prediction/" + \
                          "download_links/raptorx_download_link.txt",
         license_key_in="nodes/utils/protein_structure_prediction/" + \
                        "download_links/modeller_license_key.txt"
    output:
         fasta_out=f"data/{DATASET}/annotated_pdbs_seqs.pdb",
         classes_out=f"data/{DATASET}/annotated_pdbs_classes.pdb",
         pdbs_out=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1])
    params:
         subworkflow="protein_structure_prediction",
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
         dryrun=config["dryrun"]
    resources:
         cores=-1
    script:
         "utils/subworkflow.py"

# rule util_convert_to_sdf_and_minimize:
#         input:
#              fasta_in=,
#              classes_in=,
#              pdbs_in=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
#                             seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[:10]
#         output:
#              fasta_out=,
#              classes_out,
#              pdbs_out=expand(f"data/{DATASET}/pdb/{{seq_name}}_minimized.pdb",
#                             seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[:10]
#         params:
#              subworkflow="convert_to_sdf_and_minimize",
#              snakefile="nodes/utils/convert_to_sdf_and_minimize/Snakefile",
#              configfile="nodes/utils/convert_to_sdf_and_minimize/config.yaml",
#              dryrun=config["dryrun"]
#         resources:
#              cores=1  # can be omitted, uses one core per default
#         script:
#              "utils/subworkflow.py"

rule encoding_qsar:
    input:
         fasta_in=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdb_dir=f"data/{DATASET}/pdb/"
    output:
         csv_out=f"data/{DATASET}/csv/qsar.csv"
    params:
         subworkflow="qsar",
         snakefile="nodes/encodings/qsar/Snakefile",
         configfile="nodes/encodings/qsar/config.yaml",
         dryrun=config["dryrun"]
    resources:
         cores=1  # can be omitted, uses one core per default
    script:
         "utils/subworkflow.py"