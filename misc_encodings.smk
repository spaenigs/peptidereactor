from utils.snakemake_config import WorkflowExecuter
import pandas as pd

DATASET = config["dataset"]
CORES = 8

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

rule encoding_pssm:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=f"data/{DATASET}/csv/pssm.csv"
    params:
         snakefile="nodes/encodings/pssm/Snakefile",
         configfile="nodes/encodings/pssm/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_zscale:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/zscale.csv"
    params:
         snakefile="nodes/encodings/zscale/Snakefile",
         configfile="nodes/encodings/zscale/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_tpc:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/tpc.csv"
    params:
         snakefile="nodes/encodings/tpc/Snakefile",
         configfile="nodes/encodings/tpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_gtpc:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/gtpc.csv"
    params:
         snakefile="nodes/encodings/gtpc/Snakefile",
         configfile="nodes/encodings/gtpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_gdpc:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/gdpc.csv"
    params:
         snakefile="nodes/encodings/gdpc/Snakefile",
         configfile="nodes/encodings/gdpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_gaac:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/gaac.csv"
    params:
         snakefile="nodes/encodings/gaac/Snakefile",
         configfile="nodes/encodings/gaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_egaac:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/egaac/egaac_window_{{window_val}}.csv",
                        window_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/egaac/Snakefile",
         configfile="nodes/encodings/egaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_dpc:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/dpc.csv"
    params:
         snakefile="nodes/encodings/dpc/Snakefile",
         configfile="nodes/encodings/dpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_dde:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/dde.csv"
    params:
         snakefile="nodes/encodings/dde/Snakefile",
         configfile="nodes/encodings/dde/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ctdt:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/ctdt.csv"
    params:
         snakefile="nodes/encodings/ctdt/Snakefile",
         configfile="nodes/encodings/ctdt/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ctdd:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/ctdd.csv"
    params:
         snakefile="nodes/encodings/ctdd/Snakefile",
         configfile="nodes/encodings/ctdd/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ctdc:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/ctdc.csv"
    params:
         snakefile="nodes/encodings/ctdc/Snakefile",
         configfile="nodes/encodings/ctdc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_blosum62:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/blosum62.csv"
    params:
         snakefile="nodes/encodings/blosum62/Snakefile",
         configfile="nodes/encodings/blosum62/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_binary:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs_msa.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/binary.csv"
    params:
         snakefile="nodes/encodings/binary/Snakefile",
         configfile="nodes/encodings/binary/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_aaindex:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/aaindex/aaindex_{{aaindex}}.csv",
                         aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/aaindex/Snakefile",
         configfile="nodes/encodings/aaindex/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")



rule encoding_aac:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/aac.csv"
    params:
         snakefile="nodes/encodings/aac/Snakefile",
         configfile="nodes/encodings/aac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_fft:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         csv_in=expand(f"data/{DATASET}/csv/aaindex/aaindex_{{aaindex}}.csv",
                       aaindex=get_aaindex())
    output:
         csv_out=expand(f"data/{DATASET}/csv/fft/fft_{{aaindex}}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/fft/Snakefile",
         configfile="nodes/encodings/distance_distribution/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores 8 \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_cgr:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/cgr/cgr_res_{{resolution}}_sf_{{sfactor}}.csv",
                        resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713])
    params:
         snakefile="nodes/encodings/cgr/Snakefile",
         configfile="nodes/encodings/cgr/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ctriad:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/ctriad.csv"
    params:
         snakefile="nodes/encodings/ctriad/Snakefile",
         configfile="nodes/encodings/ctriad/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")