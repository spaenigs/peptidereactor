from utils.snakemake_config import WorkflowExecuter

DATASET = config["dataset"]
CORES = 8

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
             shell(f"""snakemake -s {{params.snakefile}} {{output.length_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

