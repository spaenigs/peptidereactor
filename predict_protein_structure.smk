import os
from modlamp.core import read_fasta

config["global_workdir"] = os.getcwd() + "/"

DATASET = config["dataset"]

rule all:
    input:
         # expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
         #        seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
         #        ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"])
         # sdfs_out=expand(f"data/{DATASET}/sdf/{{seq_name}}.sdf",
         #                 seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[4]
         expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[-1]

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
         pdbs_out=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1])[-1]
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
#              pdbs_in=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
#                             seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[4]
#         output:
#              sdfs_out=expand(f"data/{DATASET}/sdf/{{seq_name}}.sdf",
#                             seq_name=read_fasta(f"data/{DATASET}/seqs.fasta")[1])[4]
#         params:
#              subworkflow="convert_to_sdf_and_minimize",
#              snakefile="nodes/utils/convert_to_sdf_and_minimize/Snakefile",
#              configfile="nodes/utils/convert_to_sdf_and_minimize/config.yaml",
#              dryrun=config["dryrun"]
#         resources:
#              cores=-1  # can be omitted, uses one core per default
#         script:
#              "utils/subworkflow.py"