from proteinreactor.workflow_executer import WorkflowExecuter

CORES = config["cores"]

rule all:
    input:
         "data/bachem/plots/psekraac_filtered/",
         "data/bachem/csv/non_empty/filtered/psekraac/",
         "data/bachem/csv/non_empty/filtered/aaindex/",
         "data/bachem/csv/non_empty/filtered/rest/",
         "data/bachem/csv/non_empty/filtered/total/"

rule filter_psekraac:
    input:
         base_dir_in="data/bachem/csv/non_empty/all/"
    output:
         plot_dir_out=directory("data/bachem/plots/psekraac_filtered/"),
         csv_dir_out=directory("data/bachem/csv/non_empty/filtered/psekraac/")
    params:
         snakefile="nodes/filter/psekraac/Snakefile",
         configfile="nodes/filter/psekraac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output),
                               params.configfile,
                               cores=CORES,
                               window_lengths=[8,11,13,15,17,20]) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule filter_aaindex:
    input:
         base_dir_in="data/bachem/csv/non_empty/all/"
    output:
         csv_dir_out=directory("data/bachem/csv/non_empty/filtered/aaindex/")
    params:
         snakefile="nodes/filter/aaindex/Snakefile",
         configfile="nodes/filter/aaindex/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# https://academic.oup.com/bioinformatics/article/21/8/1509/249540#2950870
rule filter_curse_of_dimensionality:
    input:
         "data/bachem/csv/non_empty/all/",
         "data/bachem/csv/non_empty/filtered/psekraac/",
         "data/bachem/csv/non_empty/filtered/aaindex/"
    output:
         directory("data/bachem/csv/non_empty/filtered/rest/")
    run:
         from glob import glob
         import pandas as pd
         enc_types = ["aaindex", "fft", "waac", "flgc", "fldpc", "psekraac"]
         for p in glob(str(input[0]) + "*.csv"):
             if sum([1 if enc_type in p else 0 for enc_type in enc_types]) > 0:
                 continue
             else:
                 df = pd.read_csv(p, index_col=0, engine="c")
                 nrows, ncols = df.shape
                 if ncols - 1 <= nrows:
                     shell(f"cp {p} {str(output)}")

rule collect:
    input:
         "data/bachem/csv/non_empty/filtered/psekraac/",
         "data/bachem/csv/non_empty/filtered/aaindex/",
         "data/bachem/csv/non_empty/filtered/rest/"
    output:
         directory("data/bachem/csv/non_empty/filtered/total/")
    run:
         from glob import glob

         def copy_files(file_list):
             for f in file_list:
                 shell(f"cp {f} {str(output)}")

         copy_files(
             glob(str(input[0]) + "*.csv") +
             glob(str(input[1]) + "*.csv") +
             glob(str(input[2]) + "*.csv")
         )







