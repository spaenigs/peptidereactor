import secrets


def _get_header(token):
    return f'''
rule benchmark_cross_validation_single_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(csv_in, csv_dir_out):
    return f'''\
    input:
        csv_in={csv_in}
    output:
        cv_dir_out=directory("{csv_dir_out}")
    threads:
         1000
    params:
        snakefile="nodes/benchmark/cross_validation/single/Snakefile",
        configfile="nodes/benchmark/cross_validation/single/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(csv_in, csv_dir_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_cross_validation_single_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_in, csv_dir_out)
    return rule
