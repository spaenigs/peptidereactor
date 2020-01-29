import os, sys

sys.path.append(os.getcwd())

from proteinreactor.workflow_executer import WorkflowExecuter

DATASET = "bachem"
DATASETS = ["bachem_window_length_8"]
DATASETS_COMPLETE = [n + "_complete" for n in DATASETS]
CORES = int(config["cores"])
TOKEN = config["token"]

rule all:
    input:
         config["fasta_anno_out"],
         config["fasta_anno_msa_out"],
         config["classes_anno_idx_out"],
         config["profile_dir"],
         config["fasta_anno_pdbs_out"],
         config["classes_anno_pdbs_idx_out"],
         config["pdb_out"],
         config["asa_out"],
         config["ta_out"],
         config["ssec_out"],
         config["sseb_out"],
         config["disorder_out"],
         config["disorderb_out"],
         config["disorderc_out"],
         config["qsar_out"],
         config["electrostatic_hull_out"],
         config["distance_distribution_out"],
         config["delaunay_out"]

rule utils_validate_sequence_names:
    input:
         fasta_in=config["fasta_in"]
    output:
         fasta_out=temp(f"data/temp/{TOKEN}/seqs_validated.fasta")
    params:
         snakefile="nodes/utils/validate_sequence_names/Snakefile",
         configfile="nodes/utils/validate_sequence_names/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule util_secondary_structure_profile:
    input:
         fasta_in=f"data/temp/{TOKEN}/seqs_validated.fasta",
         fasta_msa_in=config["fasta_msa_in"],
         classes_in=config["classes_idx_in"],
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
         classes_anno=config["classes_anno_idx_out"],
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
         fasta_in=f"data/temp/{TOKEN}/seqs_validated.fasta",
         classes_in=config["classes_idx_in"],
         download_link_in=\
             "nodes/utils/protein_structure_prediction/download_links/raptorx_download_link.txt",
         license_key_in=\
             "nodes/utils/protein_structure_prediction/download_links/modeller_license_key.txt"
    output:
         fasta_out=config["fasta_anno_pdbs_out"],
         classes_out=config["classes_anno_pdbs_idx_out"],
         pdbs_out=directory(config["pdb_out"])
    priority:
         1000
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_asa_windowed:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["asa_out"]
    params:
         snakefile="nodes/encodings/asa/asa_windowed.smk",
         configfile="nodes/encodings/asa/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_ta_windowed:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["ta_out"]
    params:
         snakefile="nodes/encodings/ta/ta_windowed.smk",
         configfile="nodes/encodings/ta/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_ssec_windowed:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["ssec_out"]
    params:
         snakefile="nodes/encodings/ssec/ssec_windowed.smk",
         configfile="nodes/encodings/ssec/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_sseb_windowed:
    input:
         fasta_in=config["fasta_anno_msa_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["sseb_out"]
    params:
         snakefile="nodes/encodings/sseb/sseb_windowed.smk",
         configfile="nodes/encodings/sseb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_disorder_windowed:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorder_out"]
    params:
         snakefile="nodes/encodings/disorder/disorder_windowed.smk",
         configfile="nodes/encodings/disorder/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_disorderb_windowed:
    input:
         fasta_in=config["fasta_anno_msa_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorderb_out"]
    params:
         snakefile="nodes/encodings/disorderb/disorderb_windowed.smk",
         configfile="nodes/encodings/disorderb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_disorderc_windowed:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_idx_in=config["classes_anno_idx_out"],
         classes_in=config["classes_in"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorderc_out"]
    params:
         snakefile="nodes/encodings/disorderc/disorderc_windowed.smk",
         configfile="nodes/encodings/disorderc/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_qsar_windowed:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_idx_in=config["classes_anno_pdbs_idx_out"],
         classes_in=config["classes_in"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["qsar_out"]
    params:
         snakefile="nodes/encodings/qsar/qsar_windowed.smk",
         configfile="nodes/encodings/qsar/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_electrostatic_hull_windowed:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_idx_in=config["classes_anno_pdbs_idx_out"],
         classes_in=config["classes_in"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["electrostatic_hull_out"]
    params:
         snakefile="nodes/encodings/electrostatic_hull/electrostatic_hull_windowed.smk",
         configfile="nodes/encodings/electrostatic_hull/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_distance_distribution_windowed:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_idx_in=config["classes_anno_pdbs_idx_out"],
         classes_in=config["classes_in"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["distance_distribution_out"]
    params:
         snakefile="nodes/encodings/distance_distribution/distance_distribution_windowed.smk",
         configfile="nodes/encodings/distance_distribution/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule encoding_delaunay_windowed:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_idx_in=config["classes_anno_pdbs_idx_out"],
         classes_in=config["classes_in"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["delaunay_out"]
    params:
         snakefile="nodes/encodings/delaunay/delaunay_windowed.smk",
         configfile="nodes/encodings/delaunay/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")


