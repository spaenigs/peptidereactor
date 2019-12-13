from utils.snakemake_config import WorkflowExecuter

DATASET = config["dataset"]
CORES = 8

rule util_window_length:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta"
    output:
         length_out=expand(f"data/{DATASET}/misc/{{encoding}}.yaml",
                           encoding=["ksctriad", "moran", "nmbroto", "geary",
                                     "qsorder", "socnumber", "eaac", "cksaagp",
                                     "cksaap", "apaac", "paac"])
    params:
         snakefile="nodes/utils/window_length/Snakefile",
         configfile="nodes/utils/window_length/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.length_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")