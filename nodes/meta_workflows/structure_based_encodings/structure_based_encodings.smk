import os, sys

sys.path.append(os.getcwd())

from peptidereactor.workflow_executer import WorkflowExecuter

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

rule encoding_pssm:
    input:
         fasta_in=config["fasta_anno"],
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
         fasta_in=config["fasta_anno"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["asa_out"]
    params:
         snakefile="nodes/encodings/asa/asa.smk",
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
         fasta_in=config["fasta_anno"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["ta_out"]
    params:
         snakefile="nodes/encodings/ta/ta.smk",
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
         fasta_in=config["fasta_anno"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["ssec_out"]
    params:
         snakefile="nodes/encodings/ssec/ssec.smk",
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
         fasta_in=config["fasta_anno_msa"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["sseb_out"]
    params:
         snakefile="nodes/encodings/sseb/sseb.smk",
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
         fasta_in=config["fasta_anno"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorder_out"]
    params:
         snakefile="nodes/encodings/disorder/disorder.smk",
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
         fasta_in=config["fasta_anno_msa"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorderb_out"]
    params:
         snakefile="nodes/encodings/disorderb/disorderb.smk",
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
         fasta_in=config["fasta_anno"],
         classes_in=config["classes_anno"],
         profile=config["profile_dir"]
    output:
         csv_out=config["disorderc_out"]
    params:
         snakefile="nodes/encodings/disorderc/disorderc.smk",
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
         fasta_in=config["fasta_anno_pdbs"],
         classes_in=config["classes_anno_pdbs"],
         pdb_dir=config["pdb_dir"]
    output:
         csv_out=config["qsar_out"]
    params:
         snakefile="nodes/encodings/qsar/qsar.smk",
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
         fasta_in=config["fasta_anno_pdbs"],
         classes_in=config["classes_anno_pdbs"],
         pdb_dir=config["pdb_dir"]
    output:
         csv_out=config["electrostatic_hull_out"]
    params:
         snakefile="nodes/encodings/electrostatic_hull/electrostatic_hull.smk",
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
         fasta_in=config["fasta_anno_pdbs"],
         classes_in=config["classes_anno_pdbs"],
         pdb_dir=config["pdb_dir"]
    output:
         csv_out=config["distance_distribution_out"]
    params:
         snakefile="nodes/encodings/distance_distribution/distance_distribution.smk",
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
         fasta_in=config["fasta_anno_pdbs"],
         classes_in=config["classes_anno_pdbs"],
         pdb_dir=config["pdb_dir"]
    output:
         csv_out=config["delaunay_out"]
    params:
         snakefile="nodes/encodings/delaunay/delaunay.smk",
         configfile="nodes/encodings/delaunay/config.yaml"
    run:
         check_empty(path_to_fasta=input.fasta_in,
                     path_to_csv_out=output.csv_out,
                     dict_input=dict(input),
                     dict_output=dict(output),
                     params_snakefile=params.snakefile,
                     params_configfile= params.configfile)
