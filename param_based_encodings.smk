from utils.snakemake_config import WorkflowExecuter
import sys
import yaml

DATASET = config["dataset"]
CORES = 8

max_vals = {}
for encoding in ["ksctriad", "moran", "nmbroto", "geary","qsorder",
                 "socnumber", "eaac", "cksaagp", "cksaap", "apaac", "paac"]:
    try:
        with open(f"data/{DATASET}/misc/{encoding}.yaml") as f:
            max_vals[encoding] = yaml.safe_load(f) + 1  # range is exclusive
    except FileNotFoundError:
        sys.exit("""
        Please run node window_length beforehand (or set them manually).
        See, e.g., apps/iFeature/codes/KSCTriad.py for details.
        """)

rule encoding_cksaagp:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/cksaagp.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/cksaagp/cksaagp_gap_{{gap_val}}.csv",
                        gap_val=list(range(1, max_vals["cksaagp"])))
    params:
         snakefile="nodes/encodings/cksaagp/Snakefile",
         configfile="nodes/encodings/cksaagp/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_socnumber:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/socnumber.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/socnumber/socnumber_nlag_{{nlag_val}}.csv",
                        nlag_val=list(range(1, max_vals["socnumber"])))
    params:
         snakefile="nodes/encodings/socnumber/Snakefile",
         configfile="nodes/encodings/socnumber/config.yaml"
    resources:
         cores=4
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_qsorder:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/qsorder.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/qsorder/qsorder_nlag_{{nlag_val}}.csv",
                        nlag_val=list(range(1, max_vals["qsorder"])))
    params:
         snakefile="nodes/encodings/qsorder/Snakefile",
         configfile="nodes/encodings/qsorder/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_nmbroto:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/nmbroto.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/nmbroto/nmbroto_nlag_{{nlag_val}}.csv",
                        nlag_val=list(range(1, max_vals["nmbroto"])))
    params:
         snakefile="nodes/encodings/nmbroto/Snakefile",
         configfile="nodes/encodings/nmbroto/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_moran:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/moran.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/moran/moran_nlag_{{nlag_val}}.csv",
                        nlag_val=list(range(1, max_vals["moran"])))
    params:
         snakefile="nodes/encodings/moran/Snakefile",
         configfile="nodes/encodings/moran/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ksctriad:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/ksctriad.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ksctriad/ksctriad_gap_{{gap_val}}.csv",
                        gap_val=list(range(1, max_vals["ksctriad"])))
    params:
         snakefile="nodes/encodings/ksctriad/Snakefile",
         configfile="nodes/encodings/ksctriad/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_geary:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/geary.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/geary/geary_nlag_{{nlag_val}}.csv",
                        nlag_val=list(range(1, max_vals["geary"])))
    params:
         snakefile="nodes/encodings/geary/Snakefile",
         configfile="nodes/encodings/geary/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_eaac:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/eaac.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/eaac/eaac_window_{{window_val}}.csv",
                        window_val=list(range(1, max_vals["eaac"])))
    params:
         snakefile="nodes/encodings/eaac/Snakefile",
         configfile="nodes/encodings/eaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")



rule encoding_cksaap:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/cksaap.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/cksaap/cksaap_gap_{{gap_val}}.csv",
                        gap_val=list(range(1, max_vals["cksaap"])))
    params:
         snakefile="nodes/encodings/cksaap/Snakefile",
         configfile="nodes/encodings/cksaap/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_apaac:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/apaac.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/apaac/apaac_lambda_{{lambda_val}}.csv",
                        lambda_val=list(range(1, max_vals["apaac"])))
    params:
         snakefile="nodes/encodings/apaac/Snakefile",
         configfile="nodes/encodings/apaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_paac:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/paac.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/paac/paac_lambda_{{lambda_val}}.csv",
                        lambda_val=list(range(1, max_vals["paac"])))
    params:
         snakefile="nodes/encodings/paac/Snakefile",
         configfile="nodes/encodings/paac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")