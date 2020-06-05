import secrets


def _get_header(token):
    return f'''
rule filter_aggregate_directories_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(dirs_in, dir_out):
    return f'''\
    input:
         dirs_in={dirs_in}
    output:
         dir_out=directory("{dir_out}")
    threads:
         1000              
    params:
         snakefile="nodes/filter/aggregate_directories/Snakefile",
         configfile="nodes/filter/aggregate_directories/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(dirs_in, dir_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}filter_aggregate_directories_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(dirs_in, dir_out)
    return rule
