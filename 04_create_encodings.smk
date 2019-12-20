from utils.snakemake_config import WorkflowExecuter
import pandas as pd
import yaml
import sys

DATASET = config["dataset"]
CORES = int(config["cores"])

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

def get_max_vals(encoding):
    try:
        with open(f"data/{DATASET}/misc/{encoding}.yaml") as f:
            return yaml.safe_load(f) + 1  # range is exclusive
    except FileNotFoundError:
        sys.exit("""
        Please run node window_length beforehand (or set them manually).
        See, e.g., apps/iFeature/codes/KSCTriad.py for details.
        """)

def get_max_dim_size(ngram):
    ngram_type, size = list(ngram)
    try:
        with open(f"data/{DATASET}/misc/ngram_{ngram_type}{size}.yaml") as f:
            return yaml.safe_load(f)  # range is exclusive
    except FileNotFoundError:
        sys.exit("""
        Please run node dim_size beforehand (or set dimension manually): min(len(shape[0], shape[1]).
        """)

rule all:
    input:
        expand(f"data/{DATASET}/csv/aaindex/aaindex_{{aaindex}}.csv",
               aaindex=get_aaindex()),
        expand(f"data/{DATASET}/csv/apaac/apaac_lambda_{{lambda_val}}.csv",
               lambda_val=list(range(1, get_max_vals("apaac")))),
        expand(f"data/{DATASET}/csv/cksaagp/cksaagp_gap_{{gap_val}}.csv",
               gap_val=list(range(1, get_max_vals("cksaagp")))),
        expand(f"data/{DATASET}/csv/cksaap/cksaap_gap_{{gap_val}}.csv",
               gap_val=list(range(1, get_max_vals("cksaap")))),
        expand(f"data/{DATASET}/csv/eaac/eaac_window_{{window_val}}.csv",
               window_val=list(range(1, get_max_vals("eaac")))),
        expand(f"data/{DATASET}/csv/egaac/egaac_window_{{window_val}}.csv",
               window_val=list(range(1, 31))),
        expand(f"data/{DATASET}/csv/geary/geary_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("geary")))),
        expand(f"data/{DATASET}/csv/ksctriad/ksctriad_gap_{{gap_val}}.csv",
               gap_val=list(range(1, get_max_vals("ksctriad")))),
        expand(f"data/{DATASET}/csv/moran/moran_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("moran")))),
        expand(f"data/{DATASET}/csv/nmbroto/nmbroto_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("nmbroto")))),
        expand(f"data/{DATASET}/csv/paac/paac_lambda_{{lambda_val}}.csv",
               lambda_val=list(range(1, get_max_vals("paac")))),
        expand(f"data/{DATASET}/csv/qsorder/qsorder_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("qsorder")))),
        expand(f"data/{DATASET}/csv/socnumber/socnumber_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("socnumber")))),

        expand(f"data/{DATASET}/csv/psekraac_type1/psekraac_type1_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type2/psekraac_type2_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type3A/psekraac_type3A_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type3B/psekraac_type3B_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type4/psekraac_type4_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type5/psekraac_type5_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type6A/psekraac_type6A_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type6B/psekraac_type6B_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type6C/psekraac_type6C_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type7/psekraac_type7_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type8/psekraac_type8_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type9/psekraac_type9_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type10/psekraac_type10_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type11/psekraac_type11_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type12/psekraac_type12_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type13/psekraac_type13_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type14/psekraac_type14_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type15/psekraac_type15_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type16/psekraac_type16_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/fft/fft_{{aaindex}}.csv", aaindex=get_aaindex()),
        expand(f"data/{DATASET}/csv/cgr/cgr_res_{{resolution}}_sf_{{sfactor}}.csv",
               resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),

        expand(f"data/{DATASET}/csv/electrostatic_hull/electrostatic_hull_{{distance}}.csv",
               distance=[0,3,6,9,12]),
        expand(f"data/{DATASET}/csv/delaunay/delaunay_{{algorithm}}.csv",
               algorithm=["average_distance", "total_distance", "cartesian_product",
                          "number_instances", "frequency_instances"]),
        expand(f"data/{DATASET}/csv/waac/waac_{{aaindex}}.csv",
               aaindex=get_aaindex()),
        expand(f"data/{DATASET}/csv/flgc/flgc_{{aaindex}}.csv",
               aaindex=get_aaindex()),
        expand(f"data/{DATASET}/csv/fldpc/fldpc_{{aaindex}}.csv",
               aaindex=get_aaindex()),

        expand(f"data/{DATASET}/csv/ngram_a2/ngram_a2_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a2"))),
        expand(f"data/{DATASET}/csv/ngram_a2/ngram_a2_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a2"))),
        expand(f"data/{DATASET}/csv/ngram_a2/ngram_a2_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a2"))),

        expand(f"data/{DATASET}/csv/ngram_a3/ngram_a3_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a3"))),
        expand(f"data/{DATASET}/csv/ngram_a3/ngram_a3_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a3"))),
        expand(f"data/{DATASET}/csv/ngram_a3/ngram_a3_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a3"))),

        expand(f"data/{DATASET}/csv/ngram_e2/ngram_e2_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e2"))),
        expand(f"data/{DATASET}/csv/ngram_e2/ngram_e2_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e2"))),
        expand(f"data/{DATASET}/csv/ngram_e2/ngram_e2_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e2"))),

        expand(f"data/{DATASET}/csv/ngram_e3/ngram_e3_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e3"))),
        expand(f"data/{DATASET}/csv/ngram_e3/ngram_e3_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e3"))),
        expand(f"data/{DATASET}/csv/ngram_e3/ngram_e3_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e3"))),

        expand(f"data/{DATASET}/csv/ngram_s2/ngram_s2_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s2"))),
        expand(f"data/{DATASET}/csv/ngram_s2/ngram_s2_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s2"))),
        expand(f"data/{DATASET}/csv/ngram_s2/ngram_s2_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s2"))),

        expand(f"data/{DATASET}/csv/ngram_s3/ngram_s3_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s3"))),
        expand(f"data/{DATASET}/csv/ngram_s3/ngram_s3_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s3"))),
        expand(f"data/{DATASET}/csv/ngram_s3/ngram_s3_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s3"))),

        expand(f"data/{DATASET}/csv/distance_frequency/distance_frequency_dn_{{nterminal}}_dc_{{cterminal}}.csv",
               nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),

        f"data/{DATASET}/csv/disorder.csv",
        f"data/{DATASET}/csv/disorderb.csv",
        f"data/{DATASET}/csv/disorderc.csv",
        f"data/{DATASET}/csv/aac.csv",
        f"data/{DATASET}/csv/binary.csv",
        f"data/{DATASET}/csv/blosum62.csv",
        f"data/{DATASET}/csv/ctdc.csv",
        f"data/{DATASET}/csv/ctdd.csv",
        f"data/{DATASET}/csv/ctdt.csv",
        f"data/{DATASET}/csv/ctriad.csv",
        f"data/{DATASET}/csv/dde.csv",
        f"data/{DATASET}/csv/dpc.csv",
        f"data/{DATASET}/csv/gaac.csv",
        f"data/{DATASET}/csv/gdpc.csv",
        f"data/{DATASET}/csv/gtpc.csv",
        f"data/{DATASET}/csv/tpc.csv",
        f"data/{DATASET}/csv/zscale.csv",
        f"data/{DATASET}/csv/pssm.csv",
        f"data/{DATASET}/csv/sseb.csv",
        f"data/{DATASET}/csv/ssec.csv",
        f"data/{DATASET}/csv/ta.csv",
        f"data/{DATASET}/csv/asa.csv",
        f"data/{DATASET}/csv/blomap.csv",
        f"data/{DATASET}/csv/distance_distribution.csv",
        f"data/{DATASET}/csv/qsar.csv"

########################################################################################################################
############################################## MISC ENCODINGS ##########################################################
########################################################################################################################

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
         configfile="nodes/encodings/fft/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_waac:
    input:
         csv_in=f"data/{DATASET}/csv/aac.csv",
         aaindex_in="apps/iFeature/data/AAindex.tsv"
    output:
         csv_out=expand(f"data/{DATASET}/csv/waac/waac_{{aaindex}}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/waac/Snakefile",
         configfile="nodes/encodings/waac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_flgc:
    input:
         csv_in=f"data/{DATASET}/csv/aac.csv",
         aaindex_in="apps/iFeature/data/AAindex.tsv"
    output:
         csv_out=expand(f"data/{DATASET}/csv/flgc/flgc_{{aaindex}}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/flgc/Snakefile",
         configfile="nodes/encodings/flgc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_fldpc:
    input:
         csv_in=f"data/{DATASET}/csv/dpc.csv",
         aaindex_in="apps/iFeature/data/AAindex.tsv"
    output:
         csv_out=expand(f"data/{DATASET}/csv/fldpc/fldpc_{{aaindex}}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/fldpc/Snakefile",
         configfile="nodes/encodings/fldpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ngram_a2:
    input:
         csv_in=f"data/{DATASET}/csv/dpc.csv",
         length_in=f"data/{DATASET}/misc/ngram_a2.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ngram_a2/ngram_a2_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("a2"))),
         lsv_out=expand(f"data/{DATASET}/csv/ngram_a2/ngram_a2_lsv_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("a2"))),
         sv_out=expand(f"data/{DATASET}/csv/ngram_a2/ngram_a2_sv_{{dim}}.csv",
                       dim=range(1, get_max_dim_size("a2"))),
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ngram_a3:
    input:
         csv_in=f"data/{DATASET}/csv/tpc.csv",
         length_in=f"data/{DATASET}/misc/ngram_a3.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ngram_a3/ngram_a3_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("a3"))),
         lsv_out=expand(f"data/{DATASET}/csv/ngram_a3/ngram_a3_lsv_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("a3"))),
         sv_out=expand(f"data/{DATASET}/csv/ngram_a3/ngram_a3_sv_{{dim}}.csv",
                       dim=range(1, get_max_dim_size("a3"))),
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ngram_e2:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/ngram_e2.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ngram_e2/ngram_e2_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("e2"))),
         lsv_out=expand(f"data/{DATASET}/csv/ngram_e2/ngram_e2_lsv_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("e2"))),
         sv_out=expand(f"data/{DATASET}/csv/ngram_e2/ngram_e2_sv_{{dim}}.csv",
                       dim=range(1, get_max_dim_size("e2")))
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ngram_e3:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/ngram_e3.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ngram_e3/ngram_e3_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("e3"))),
         lsv_out=expand(f"data/{DATASET}/csv/ngram_e3/ngram_e3_lsv_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("e3"))),
         sv_out=expand(f"data/{DATASET}/csv/ngram_e3/ngram_e3_sv_{{dim}}.csv",
                       dim=range(1, get_max_dim_size("e3")))
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ngram_s2:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/ngram_s2.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ngram_s2/ngram_s2_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("s2"))),
         lsv_out=expand(f"data/{DATASET}/csv/ngram_s2/ngram_s2_lsv_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("s2"))),
         sv_out=expand(f"data/{DATASET}/csv/ngram_s2/ngram_s2_sv_{{dim}}.csv",
                       dim=range(1, get_max_dim_size("s2")))
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_ngram_s3:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/ngram_s3.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/ngram_s3/ngram_s3_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("s3"))),
         lsv_out=expand(f"data/{DATASET}/csv/ngram_s3/ngram_s3_lsv_{{dim}}.csv",
                        dim=range(1, get_max_dim_size("s3"))),
         sv_out=expand(f"data/{DATASET}/csv/ngram_s3/ngram_s3_sv_{{dim}}.csv",
                       dim=range(1, get_max_dim_size("s3")))
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
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

rule encoding_distance_frequency:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/distance_frequency/" + \
                        f"distance_frequency_dn_{{nterminal}}_dc_{{cterminal}}.csv",
                        nterminal=[5, 10, 20, 50, 100],
                        cterminal=[5, 10, 20, 50, 100])
    params:
         snakefile="nodes/encodings/distance_frequency/Snakefile",
         configfile="nodes/encodings/distance_frequency/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_blomap:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=f"data/{DATASET}/csv/blomap.csv"
    params:
         snakefile="nodes/encodings/blomap/Snakefile",
         configfile="nodes/encodings/blomap/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

########################################################################################################################
########################################### PARAM_BASED ENCODINGS ######################################################
########################################################################################################################

rule encoding_cksaagp:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt",
         length_in=f"data/{DATASET}/misc/cksaagp.yaml"
    output:
         csv_out=expand(f"data/{DATASET}/csv/cksaagp/cksaagp_gap_{{gap_val}}.csv",
                        gap_val=list(range(1, get_max_vals("cksaagp"))))
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
                        nlag_val=list(range(1, get_max_vals("socnumber"))))
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
                        nlag_val=list(range(1, get_max_vals("qsorder"))))
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
                        nlag_val=list(range(1, get_max_vals("nmbroto"))))
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
                        nlag_val=list(range(1, get_max_vals("moran"))))
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
                        gap_val=list(range(1, get_max_vals("ksctriad"))))
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
                        nlag_val=list(range(1, get_max_vals("geary"))))
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
                        window_val=list(range(1, get_max_vals("eaac"))))
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
                        gap_val=list(range(1, get_max_vals("cksaap"))))
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
                        lambda_val=list(range(1, get_max_vals("apaac"))))
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
                        lambda_val=list(range(1, get_max_vals("paac"))))
    params:
         snakefile="nodes/encodings/paac/Snakefile",
         configfile="nodes/encodings/paac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

########################################################################################################################
############################################## PSEKRAAC ENCODINGS ######################################################
########################################################################################################################

rule encoding_psekraac_type16:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type16/psekraac_type16_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type16/Snakefile",
         configfile="nodes/encodings/psekraac_type16/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type15:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type15/psekraac_type15_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type15/Snakefile",
         configfile="nodes/encodings/psekraac_type15/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type14:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type14/psekraac_type14_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type14/Snakefile",
         configfile="nodes/encodings/psekraac_type14/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type13:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type13/psekraac_type13_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type13/Snakefile",
         configfile="nodes/encodings/psekraac_type13/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type12:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type12/psekraac_type12_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type12/Snakefile",
         configfile="nodes/encodings/psekraac_type12/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type11:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type11/psekraac_type11_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type11/Snakefile",
         configfile="nodes/encodings/psekraac_type11/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type10:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type10/psekraac_type10_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type10/Snakefile",
         configfile="nodes/encodings/psekraac_type10/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type9:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type9/psekraac_type9_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type9/Snakefile",
         configfile="nodes/encodings/psekraac_type9/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type8:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type8/psekraac_type8_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type8/Snakefile",
         configfile="nodes/encodings/psekraac_type8/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type7:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type7/psekraac_type7_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type7/Snakefile",
         configfile="nodes/encodings/psekraac_type7/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type6C:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type6C/psekraac_type6C_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type6C/Snakefile",
         configfile="nodes/encodings/psekraac_type6C/config.yaml"
    resources:
         cores=8
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type6B:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type6B/psekraac_type6B_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type6B/Snakefile",
         configfile="nodes/encodings/psekraac_type6B/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type6A:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type6A/psekraac_type6A_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type6A/Snakefile",
         configfile="nodes/encodings/psekraac_type6A/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type5:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type5/psekraac_type5_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type5/Snakefile",
         configfile="nodes/encodings/psekraac_type5/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type4:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type4/psekraac_type4_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type4/Snakefile",
         configfile="nodes/encodings/psekraac_type4/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type3B:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type3B/psekraac_type3B_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type3B/Snakefile",
         configfile="nodes/encodings/psekraac_type3B/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type3A:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type3A/psekraac_type3A_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type3A/Snakefile",
         configfile="nodes/encodings/psekraac_type3A/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type2:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type2/psekraac_type2_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type2/Snakefile",
         configfile="nodes/encodings/psekraac_type2/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule encoding_psekraac_type1:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta",
         classes_in=f"data/{DATASET}/annotated_classes.txt"
    output:
         csv_out=expand(f"data/{DATASET}/csv/psekraac_type1/psekraac_type1_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type1/Snakefile",
         configfile="nodes/encodings/psekraac_type1/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.csv_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

########################################################################################################################
############################################ STRUCTURE-BASED ENCODINGS #################################################
########################################################################################################################

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