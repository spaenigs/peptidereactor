from utils.snakemake_config import WorkflowExecuter

DATASET = config["dataset"]
CORES = 8

rule encoding_asa:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/asa.csv"
    params:
         snakefile="nodes/encodings/asa/Snakefile",
         configfile="nodes/encodings/asa/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ta:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/ta.csv"
    params:
         snakefile="nodes/encodings/ta/Snakefile",
         configfile="nodes/encodings/ta/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ssec:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/ssec.csv"
    params:
         snakefile="nodes/encodings/ssec/Snakefile",
         configfile="nodes/encodings/ssec/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_sseb:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs_msa.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/sseb.csv"
    params:
         snakefile="nodes/encodings/sseb/Snakefile",
         configfile="nodes/encodings/sseb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_disorder:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/disorder.csv"
    params:
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_disorderb:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs_msa.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/disorderb.csv"
    params:
         snakefile="nodes/encodings/disorderb/Snakefile",
         configfile="nodes/encodings/disorderb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_disorderc:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/disorderc.csv"
    params:
         snakefile="nodes/encodings/disorderc/Snakefile",
         configfile="nodes/encodings/disorderc/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_qsar:
    input:
         fasta_in=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdb_dir=f"data/{DATASET}/pdb/"
    output:
         csv_out=f"data/{DATASET}/csv/qsar.csv"
    params:
         snakefile="nodes/encodings/qsar/Snakefile",
         configfile="nodes/encodings/qsar/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_electrostatic_hull:
    input:
         fasta_in=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdb_dir=f"data/{DATASET}/pdb/"
    output:
         csv_out=expand(f"data/{DATASET}/csv/electrostatic_hull/electrostatic_hull_{{distance}}.csv",
                        distance=[0,3,6,9,12])
    params:
         snakefile="nodes/encodings/electrostatic_hull/Snakefile",
         configfile="nodes/encodings/electrostatic_hull/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_distance_distribution:
    input:
         fasta_in=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdb_dir=f"data/{DATASET}/pdb/"
    output:
         csv_out=f"data/{DATASET}/csv/distance_distribution.csv"
    params:
         snakefile="nodes/encodings/distance_distribution/Snakefile",
         configfile="nodes/encodings/distance_distribution/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_delaunay:
    input:
         fasta_in=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdb_dir=f"data/{DATASET}/pdb/"
    output:
         csv_out=expand(f"data/{DATASET}/csv/delaunay/delaunay_{{algorithm}}.csv",
                        algorithm=["average_distance", "total_distance", "cartesian_product",
                                   "number_instances", "frequency_instances"])
    params:
         snakefile="nodes/encodings/delaunay/Snakefile",
         configfile="nodes/encodings/delaunay/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")