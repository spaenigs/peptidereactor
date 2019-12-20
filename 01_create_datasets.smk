from utils.snakemake_config import WorkflowExecuter

DATASET = config["dataset"]
CORES = int(config["cores"])

rule all:
    input:
         expand("data/{dataset}/sequence_length_distribution.svg", dataset=config["dataset"]),
         f"data/{DATASET}_ds1/seqs.fasta",
         f"data/{DATASET}_ds1/classes.txt",
         f"data/{DATASET}_ds2/seqs.fasta",
         f"data/{DATASET}_ds2/classes.txt"

rule plot_sequence_length_distribution:
        input:
             fasta_in="data/{dataset}/seqs.fasta",
             classes_in="data/{dataset}/classes.txt",
             fasta_out_1="data/{dataset}_ds1/seqs.fasta",
             classes_out_1="data/{dataset}_ds1/classes.txt",
             fasta_out_2="data/{dataset}_ds2/seqs.fasta",
             classes_out_2="data/{dataset}_ds2/classes.txt"
        output:
             svg_out="data/{dataset}/sequence_length_distribution.svg"
        params:
             snakefile="nodes/plots/sequence_length_distribution/Snakefile",
             configfile="nodes/plots/sequence_length_distribution/config.yaml"
        run:
             with WorkflowExecuter(dict(input), dict(output), params.configfile):
                 shell(f"""snakemake -s {{params.snakefile}} {{output.svg_out}} \
                                --cores {CORES} \
                                --directory $PWD \
                                --configfile {{params.configfile}}""")

rule util_split_normalize:
    input:
         fasta_in="data/{dataset}/seqs.fasta",
         classes_in="data/{dataset}/classes.txt"
    output:
         fasta_out_1="data/{dataset}_ds1/seqs.fasta",
         classes_out_1="data/{dataset}_ds1/classes.txt",
         fasta_out_2="data/{dataset}_ds2/seqs.fasta",
         classes_out_2="data/{dataset}_ds2/classes.txt"
    params:
         snakefile="nodes/utils/split_normalize/Snakefile",
         configfile="nodes/utils/split_normalize/config.yaml"
    run:
         import os
         print(os.getcwd())
         print(os.listdir(os.getcwd()))
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.fasta_out_1}} {{output.classes_out_2}} {{output.fasta_out_2}} {{output.classes_out_2}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")