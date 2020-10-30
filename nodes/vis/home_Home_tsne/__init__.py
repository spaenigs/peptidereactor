import secrets


def _get_header(token):
    return f'''
rule vis_home_Home_tsne_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fastas_in, classes_in, readmes_in, benchmark_csvs_in, html_dir_out):
    return f'''
    input:
         fastas_in={fastas_in},
         classes_in={classes_in},
         readmes_in={readmes_in},
         benchmark_csvs_in={benchmark_csvs_in},    
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/home_Home_tsne/Snakefile",
         configfile="nodes/vis/home_Home_tsne/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fastas_in, classes_in, readmes_in, benchmark_csvs_in, html_dir_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_home_Home_tsne_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fastas_in, classes_in, readmes_in, benchmark_csvs_in, html_dir_out)
    return rule
