import os, sys

sys.path.append(os.getcwd())

from proteinreactor.workflow_executer import WorkflowExecuter

TOKEN = config["token"]
CORES = config["cores"]

def check_empty(path_to_fasta, path_to_csv_out,
                dict_input, dict_output,
                params_configfile, params_snakefile):
    if os.path.getsize(path_to_fasta) == 0:
        if type(path_to_csv_out) == list:
            for p in path_to_csv_out:
                shell(f"touch {p}")
        else:
            shell(f"touch {path_to_csv_out}")
    else:
        with WorkflowExecuter(dict_input, dict_output, params_configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {params_snakefile} --configfile {params_configfile}""")

rule all:
    input:
         config["fasta_anno_out"],
         config["fasta_anno_msa_out"],
         config["classes_anno"],
         config["profile_dir"],
         config["fasta_anno_pdbs_out"],
         config["classes_anno_pdbs_out"],
         config["pdb_out"],
         config["pssm_out"],
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

rule encoding_pssm:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["pssm_out"]
    params:
         snakefile="nodes/encodings/pssm/Snakefile",
         configfile="nodes/encodings/pssm/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_asa:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["asa_out"]
    params:
         snakefile="nodes/encodings/asa/Snakefile",
         configfile="nodes/encodings/asa/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_ta:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["ta_out"]
    params:
         snakefile="nodes/encodings/ta/Snakefile",
         configfile="nodes/encodings/ta/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_ssec:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["ssec_out"]
    params:
         snakefile="nodes/encodings/ssec/Snakefile",
         configfile="nodes/encodings/ssec/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_sseb:
    input:
         fasta_in=config["fasta_anno_msa_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["sseb_out"]
    params:
         snakefile="nodes/encodings/sseb/Snakefile",
         configfile="nodes/encodings/sseb/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_disorder:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorder_out"]
    params:
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_disorderb:
    input:
         fasta_in=config["fasta_anno_msa_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorderb_out"]
    params:
         snakefile="nodes/encodings/disorderb/Snakefile",
         configfile="nodes/encodings/disorderb/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_disorderc:
    input:
         fasta_in=config["fasta_anno_out"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorderc_out"]
    params:
         snakefile="nodes/encodings/disorderc/Snakefile",
         configfile="nodes/encodings/disorderc/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_qsar:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_in=config["classes_anno_pdbs_out"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["qsar_out"]
    params:
         snakefile="nodes/encodings/qsar/Snakefile",
         configfile="nodes/encodings/qsar/config.yaml",
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_electrostatic_hull:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_in=config["classes_anno_pdbs_out"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["electrostatic_hull_out"]
    params:
         snakefile="nodes/encodings/electrostatic_hull/Snakefile",
         configfile="nodes/encodings/electrostatic_hull/config.yaml"
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_distance_distribution:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_in=config["classes_anno_pdbs_out"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["distance_distribution_out"]
    params:
         snakefile="nodes/encodings/distance_distribution/Snakefile",
         configfile="nodes/encodings/distance_distribution/config.yaml"
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)

rule encoding_delaunay:
    input:
         fasta_in=config["fasta_anno_pdbs_out"],
         classes_in=config["classes_anno_pdbs_out"],
         pdb_dir=config["pdb_out"]
    output:
         csv_out=config["delaunay_out"]
    params:
         snakefile="nodes/encodings/delaunay/Snakefile",
         configfile="nodes/encodings/delaunay/config.yaml"
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)
