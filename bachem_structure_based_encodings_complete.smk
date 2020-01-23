from utils.snakemake_config import WorkflowExecuter

DATASET = "bachem"
DATASETS = ["bachem_window_length_8"]
DATASETS_COMPLETE = [n + "_complete" for n in DATASETS]
CORES = int(config["cores"])

rule all:
    input:
         # expand("data/{normalized_dataset}/seqs_validated.fasta", normalized_dataset=DATASETS_COMPLETE)
         # expand("data/{normalized_dataset}/annotated_seqs.fasta", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/annotated_seqs_msa.fasta", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/annotated_classes.txt", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/profile/", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/asa.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/ta.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/ssec.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/sseb.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/disorder.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/disorderb.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/disorderc.csv", normalized_dataset=DATASETS_COMPLETE),
         # expand("data/{normalized_dataset}/csv/qsar.csv", normalized_dataset=DATASETS_COMPLETE)
         expand("data/{normalized_dataset}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                normalized_dataset=DATASETS_COMPLETE, distance=[6])

rule utils_sliding_windows:
    input:
         series_in=f"data/{DATASET}/series.yaml"
    output:
         fastas_out=expand(f"data/{DATASET}_window_length_{{window_length}}/seqs.fasta", window_length=[8,15]),
         classes_out=expand(f"data/{DATASET}_window_length_{{window_length}}/classes.txt", window_length=[8,15])
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule utils_sliding_windows_complete:
    input:
         series_in=f"data/{DATASET}/series.yaml"
    output:
         fastas_out=expand(f"data/{DATASET}_window_length_{{window_length}}_complete/seqs.fasta", window_length=[8,15]),
         classes_out=expand(f"data/{DATASET}_window_length_{{window_length}}_complete/classes.yaml", window_length=[8,15]),
         classes_idx_out=expand(f"data/{DATASET}_window_length_{{window_length}}_complete/classes.txt", window_length=[8,15])
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows_complete.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule utils_validate_sequence_names:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta"
    output:
         fasta_out="data/{normalized_dataset}/seqs_validated.fasta"
    params:
         snakefile="nodes/utils/validate_sequence_names/Snakefile",
         configfile="nodes/utils/validate_sequence_names/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule util_multiple_sequence_alignment:
    input:
         fasta_in="data/{normalized_dataset}/seqs_validated.fasta",
         classes_in="data/{normalized_dataset}/classes.txt"
    output:
         fasta_out="data/{normalized_dataset}/seqs_msa.fasta"
    params:
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule util_secondary_structure_profile:
    input:
         fasta_in="data/{normalized_dataset}/seqs_validated.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
         uniprot90_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/uniprot90_download_link.txt",
         psipred_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/psipred_download_link.txt",
         spineXpublic_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/spineXpublic_download_link.txt",
         VSL2_download_link_in=\
             "nodes/utils/secondary_structure_profile/download_links/VSL2_download_link.txt"
    output:
         fasta_anno_out="data/{normalized_dataset}/annotated_seqs.fasta",
         fasta_anno_msa_out="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         classes_anno="data/{normalized_dataset}/annotated_classes.txt",
         profiles_out=directory("data/{normalized_dataset}/profile/")
    priority:
         1000
    params:
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule util_protein_structure_prediction:
    input:
         fasta_in="data/{normalized_dataset}/seqs_validated.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
         download_link_in=\
             "nodes/utils/protein_structure_prediction/download_links/raptorx_download_link.txt",
         license_key_in=\
             "nodes/utils/protein_structure_prediction/download_links/modeller_license_key.txt"
    output:
         fasta_out="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_out="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdbs_out=directory("data/{normalized_dataset}/pdb/")
    priority:
         1000
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

# rule encoding_binary:
#     input:
#          fasta_in="data/{normalized_dataset}/seqs_msa.fasta",
#          classes_in="data/{normalized_dataset}/classes.txt"
#     output:
#          csv_out="data/{normalized_dataset}/csv/binary.csv"
#     params:
#          snakefile="nodes/encodings/binary/Snakefile",
#          configfile="nodes/encodings/binary/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile):
#              shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_asa:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+}/asa.csv"
    params:
         snakefile="nodes/encodings/asa_windowed/asa.smk",
         configfile="nodes/encodings/asa_windowed/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_asa_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/asa.csv"
    params:
         snakefile="nodes/encodings/asa/asa_windowed.smk",
         configfile="nodes/encodings/asa/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ta_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/ta.csv"
    params:
         snakefile="nodes/encodings/ta/ta_windowed.smk",
         configfile="nodes/encodings/ta/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ssec_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/ssec.csv"
    params:
         snakefile="nodes/encodings/ssec/ssec_windowed.smk",
         configfile="nodes/encodings/ssec/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_sseb_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/sseb.csv"
    params:
         snakefile="nodes/encodings/sseb/sseb_windowed.smk",
         configfile="nodes/encodings/sseb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_disorder_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/disorder.csv"
    params:
         snakefile="nodes/encodings/disorder/disorder_windowed.smk",
         configfile="nodes/encodings/disorder/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_disorderb_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/disorderb.csv"
    params:
         snakefile="nodes/encodings/disorderb/disorderb_windowed.smk",
         configfile="nodes/encodings/disorderb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_disorderc_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_\w+}/csv/disorderc.csv"
    params:
         snakefile="nodes/encodings/disorderc/disorderc_windowed.smk",
         configfile="nodes/encodings/disorderc/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_qsar_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         csv_out="data/{normalized_dataset,\w+_window_length_\d+_complete}/csv/qsar.csv"
    params:
         snakefile="nodes/encodings/qsar/qsar_windowed.smk",
         configfile="nodes/encodings/qsar/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_electrostatic_hull_windowed:
    input:
         fasta_in="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_idx_in="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                     distance=[6])
    params:
         snakefile="nodes/encodings/electrostatic_hull/electrostatic_hull_windowed.smk",
         configfile="nodes/encodings/electrostatic_hull/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

