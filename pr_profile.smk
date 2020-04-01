from peptidereactor.workflow_executer import WorkflowExecuter, MetaWorkflowExecuter

CORES = config["cores"]

WINDOW_LENGTHS = [8,11,13,15,17,20]
DATASETS = \
    ["protease_window_length_8", "protease_window_length_8_complete"] + \
    expand(["pds_window_length_{window_length}", "pds_window_length_{window_length}_complete"],
           window_length=WINDOW_LENGTHS)

# WINDOW_LENGTHS = [13]#[8,11,15,20]
# DATASETS = ["protease_window_length_8", "pds_window_length_13"]

PDS_ALL = \
    [ds for ds in DATASETS if "complete" not in ds and "pds" in ds]

PDS_STRUC = \
    [ds for ds in DATASETS if "complete" not in ds and int(ds.split("_")[-1]) >= 12 and "pds" in ds]

PDS_STRUC_WINDOWED = \
    [ds for ds in DATASETS if "complete" in ds and "pds" in ds]

PROTEASE_ALL = \
    [ds for ds in DATASETS if "complete" not in ds and "protease" in ds]

PROTEASE_STRUC = \
    [ds for ds in DATASETS if "complete" not in ds and int(ds.split("_")[-1]) >= 12 and "protease" in ds]

PROTEASE_STRUC_WINDOWED = \
    [ds for ds in DATASETS if "complete" in ds and "protease" in ds]

rule all:
    input:
         # msa
         expand("data/{normalized_dataset}/seqs_msa.fasta",
                normalized_dataset=PDS_ALL + PDS_STRUC + PDS_STRUC_WINDOWED),
         # profiles for struc. based encodings
         expand("data/{normalized_dataset}/annotated_seqs.fasta",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         expand("data/{normalized_dataset}/annotated_classes.txt",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         expand("data/{normalized_dataset}/annotated_seqs_msa.fasta",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         expand("data/{normalized_dataset}/profile/",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         expand("data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         expand("data/{normalized_dataset}/annotated_pdbs_classes.txt",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         expand("data/{normalized_dataset}/pdb/",
                normalized_dataset=PDS_STRUC + PROTEASE_STRUC),
         # profiles for struc. based encodings (windowed)
         expand("data/{normalized_dataset}/seqs.fasta",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/annotated_seqs.fasta",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/annotated_classes.txt",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/annotated_seqs_msa.fasta",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/profile/",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/annotated_pdbs_classes.txt",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),
         expand("data/{normalized_dataset}/pdb/",
                normalized_dataset=PDS_STRUC_WINDOWED + PROTEASE_STRUC_WINDOWED),

rule utils_sliding_windows:
    input:
         series_in=f"data/pds/series.yaml"
    output:
         fastas_out=temp(expand(f"data/pds_window_length_{{window_length}}/seqstmp.fasta", window_length=WINDOW_LENGTHS)),
         classes_out=expand(f"data/pds_window_length_{{window_length}}/classes.txt", window_length=WINDOW_LENGTHS),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_sliding_windows_complete:
    input:
         series_in=f"data/pds/series.yaml"
    output:
         fastas_out=temp(expand(f"data/pds_window_length_{{window_length}}_complete/seqstmp.fasta", window_length=WINDOW_LENGTHS)),
         classes_out=expand(f"data/pds_window_length_{{window_length}}_complete/classes.yaml", window_length=WINDOW_LENGTHS),
         classes_idx_out=expand(f"data/pds_window_length_{{window_length}}_complete/classes.txt", window_length=WINDOW_LENGTHS),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows_complete.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_protein_dataset_creation:
    input:
         dataset_in="data/protease/impensData.txt",
    output:
         fasta_out=temp("data/protease_window_length_8/seqstmp.fasta"),
         classes_out="data/protease_window_length_8/classes.txt"
    params:
         snakefile="nodes/utils/protein_dataset_creation/protein_dataset_creation.smk",
         configfile="nodes/utils/protein_dataset_creation/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_protein_dataset_creation_complete:
    input:
         dataset_in="data/protease/impensData.txt",
         ids_file_in="data/protease/impens_ids.txt"
    output:
         fasta_out="data/protease_window_length_8_complete/seqs_orig.fasta",
         classes_out="data/protease_window_length_8_complete/classes_orig.txt",
         fasta_complete_out="data/protease_window_length_8_complete/seqs.fasta",
         classes_yaml_out="data/protease_window_length_8_complete/classes.yaml",
         classes_idx_out="data/protease_window_length_8_complete/classes.txt"
    params:
         snakefile="nodes/utils/protein_dataset_creation/protein_dataset_creation_complete.smk",
         configfile="nodes/utils/protein_dataset_creation/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_map_sequence_names:
    input:
         fasta_in="data/{normalized_dataset}/seqstmp.fasta"
    output:
         fasta_out="data/{normalized_dataset}/seqs.fasta",
         maps_out="data/{normalized_dataset}/misc/maps.yaml"
    params:
         snakefile="nodes/utils/map_sequence_names/Snakefile",
         configfile="nodes/utils/map_sequence_names/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="pds") as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule util_multiple_sequence_alignment:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         classes_in="data/{normalized_dataset}/classes.txt"
    output:
         fasta_out="data/{normalized_dataset}/seqs_msa.fasta"
    params:
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="pds") as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule meta_workflow_structure_based_profile:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt"
    output:
         fasta_anno_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_seqs.fasta",
         classes_anno="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?(1[2-9]|2\d)}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?(1[2-9]|2\d)}/pdb/")
    params:
         snakefile="nodes/meta_workflows/structure_based_profile/structure_based_profile.smk",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    priority:
         1000
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule meta_workflow_structure_based_profile_windowed:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_idx_in="data/{normalized_dataset}/classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml"
    output:
         fasta_anno_out="data/{normalized_dataset,.*?[a-z]}/annotated_seqs.fasta",
         classes_anno_idx_out="data/{normalized_dataset,.*?[a-z]}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?[a-z]}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?[a-z]}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?[a-z]}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_idx_out="data/{normalized_dataset,.*?[a-z]}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?[a-z]}/pdb/")
    params:
         snakefile="nodes/meta_workflows/structure_based_profile/structure_based_profile_windowed.smk",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    priority:
         1000
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
