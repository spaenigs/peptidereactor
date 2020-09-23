from peptidereactor.workflow_executer import WorkflowExecuter

CORES = 32

rule all:
    input:
         expand(['data/{dataset}/csv/cksaagp/cksaagp_gap_2.csv'], dataset=["temp/bachem"])

rule utils_window_length_05735fb1:
    input:
         fasta_in="data/{dataset}/seqs.fasta"
    output:
         length_out=['data/{dataset}/misc/cksaagp.yaml']
    threads:
         1000
    params:
         snakefile="nodes/utils/window_length/Snakefile",
         configfile="nodes/utils/window_length/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{e.snakemake} -s {params.snakefile} --configfile {params.configfile}""")

rule encoding_cksaagp_b3709e15:
    input:
         fasta_in="data/{dataset}/seqs.fasta",
         classes_in="data/{dataset}/classes.txt",
         length_in="data/{dataset}/misc/cksaagp.yaml"
    output:
         csv_out=['data/{dataset}/csv/cksaagp/cksaagp_gap_2.csv']
    threads:
         1000
    params:
         snakefile="nodes/encodings/cksaagp/Snakefile",
         configfile="nodes/encodings/cksaagp/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{e.snakemake} -s {params.snakefile} --configfile {params.configfile}""")
