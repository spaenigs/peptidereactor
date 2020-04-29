from peptidereactor.workflow_executer \
    import WorkflowExecuter

TOKEN = config["token"]
CORES = config["cores"]

rule all:
    input:
         config["fasta_anno_out"],
         config["classes_anno_idx_out"],
         config["fasta_anno_msa_out"],
         config["profile_dir"],
         config["fasta_anno_pdbs_out"],
         config["classes_anno_pdbs_idx_out"],
         config["pdb_out"]

# rule utils_validate_sequence_names:
#     input:
#          fasta_in=config["fasta_in"]
#     output:
#          fasta_out=temp(f"data/temp/{TOKEN}/seqs_validated.fasta")
#     params:
#          snakefile="nodes/utils/validate_sequence_names/Snakefile",
#          configfile="nodes/utils/validate_sequence_names/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule utils_validate_sequence_names_msa:
#     input:
#          fasta_in=config["fasta_msa_in"]
#     output:
#          fasta_out=temp(f"data/temp/{TOKEN}/seqs_msa_validated.fasta")
#     params:
#          snakefile="nodes/utils/validate_sequence_names/Snakefile",
#          configfile="nodes/utils/validate_sequence_names/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule util_secondary_structure_profile:
    input:
         fasta_in=f"data/temp/{TOKEN}/seqs_validated.fasta",
         fasta_msa_in=f"data/temp/{TOKEN}/seqs_msa_validated.fasta",
         classes_in=config["classes_idx_in"],
         uniprot90_download_link_in=\
             "nodes/utils/secondary_structure_prediction/download_links/uniprot90_download_link.txt",
         psipred_download_link_in=\
             "nodes/utils/secondary_structure_prediction/download_links/psipred_download_link.txt",
         spineXpublic_download_link_in=\
             "nodes/utils/secondary_structure_prediction/download_links/spineXpublic_download_link.txt",
         VSL2_download_link_in=\
             "nodes/utils/secondary_structure_prediction/download_links/VSL2_download_link.txt"
    output:
         fasta_anno_out=config["fasta_anno_out"],
         fasta_anno_msa_out=config["fasta_anno_msa_out"],
         classes_anno=config["classes_anno_idx_out"],
         profiles_out=directory(config["profile_dir"])
    priority:
         1000
    params:
         snakefile="nodes/utils/secondary_structure_prediction/Snakefile",
         configfile="nodes/utils/secondary_structure_prediction/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule util_protein_structure_prediction:
    input:
         fasta_in=f"data/temp/{TOKEN}/seqs_validated.fasta",
         classes_in=config["classes_idx_in"],
         download_link_in=\
             "nodes/utils/tertiary_structure_prediction/download_links/raptorx_download_link.txt",
         license_key_in=\
             "nodes/utils/tertiary_structure_prediction/download_links/modeller_license_key.txt"
    output:
         fasta_out=config["fasta_anno_pdbs_out"],
         classes_out=config["classes_anno_pdbs_idx_out"],
         pdbs_out=directory(config["pdb_out"])
    priority:
         1000
    params:
         snakefile="nodes/utils/tertiary_structure_prediction/Snakefile",
         configfile="nodes/utils/tertiary_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")