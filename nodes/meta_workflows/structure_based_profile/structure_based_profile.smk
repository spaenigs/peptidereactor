from proteinreactor.workflow_executer import WorkflowExecuter

CORES = config["cores"]

rule util_secondary_structure_profile:
    input:
         fasta_in=config["fasta_in"],
         fasta_msa_in=config["fasta_msa_in"],
         classes_in=config["classes_in"],
         uniprot90_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/uniprot90_download_link.txt",
         psipred_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/psipred_download_link.txt",
         spineXpublic_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/spineXpublic_download_link.txt",
         VSL2_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/VSL2_download_link.txt"
    output:
         fasta_anno_out=config["fasta_anno_out"],
         fasta_anno_msa_out=config["fasta_anno_msa_out"],
         classes_anno=config["classes_anno"],
         profiles_out=directory(config["profile_dir"])
    priority:
         1000
    params:
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule util_protein_structure_prediction:
    input:
         fasta_in=config["fasta_in"],
         classes_in=config["classes_in"],
         download_link_in=\
             "nodes/utils/protein_structure_prediction/download_links/raptorx_download_link.txt",
         license_key_in=\
             "nodes/utils/protein_structure_prediction/download_links/modeller_license_key.txt"
    output:
         fasta_out=config["fasta_anno_pdbs_out"],
         classes_out=config["classes_anno_pdbs_out"],
         pdbs_out=directory(config["pdb_out"])
    priority:
         1000
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")