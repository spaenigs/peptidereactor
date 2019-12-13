from utils.snakemake_config import WorkflowExecuter

DATASET = config["dataset"]
CORES = 8

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