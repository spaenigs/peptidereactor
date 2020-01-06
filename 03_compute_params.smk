from utils.snakemake_config import WorkflowExecuter

DATASET = config["dataset"]
CORES = int(config["cores"])

rule all:
    input:
         expand(f"data/{DATASET}/misc/{{encoding}}.yaml",
                encoding=["ksctriad", "moran", "nmbroto", "geary",
                          "qsorder", "socnumber", "eaac", "cksaagp",
                          "cksaap", "apaac", "paac"]),
         expand(f"data/{DATASET}/misc/ngram_{{type}}{{size}}.yaml",
                type=["a","e","s"], size=[2,3])

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
             shell(f"""snakemake -s {{params.snakefile}} \
                                 --cores {CORES} \
                                 --directory $PWD \
                                 --configfile {{params.configfile}}""")

rule util_dim_size:
    input:
         fasta_in=f"data/{DATASET}/annotated_seqs.fasta"
    output:
         length_out=expand(f"data/{DATASET}/misc/ngram_{{type}}{{size}}.yaml",
                           type=["a","e","s"], size=[2,3])
    params:
         snakefile="nodes/utils/dim_size/Snakefile",
         configfile="nodes/utils/dim_size/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} \
                                 --cores {CORES} \
                                 --directory $PWD \
                                 --configfile {{params.configfile}}""")

