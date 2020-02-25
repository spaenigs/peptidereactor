from proteinreactor.workflow_executer import WorkflowExecuter

CORES = config["cores"]

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
