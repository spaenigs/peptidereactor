import secrets


def _get_header(token):
    return f'''
rule benchmark_critical_difference_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(metrics_dir_in, cd_dir_out):
    return f'''\
    input:
        metrics_dir_in="{metrics_dir_in}"
    output:
        cd_dir_out=directory("{cd_dir_out}")
    threads:
         1000
    params:
        snakefile="nodes/benchmark/critical_difference/Snakefile",
        configfile="nodes/benchmark/critical_difference/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metrics_dir_in, cd_dir_out, benchmark_dir=None):
    """
    Computes the critical difference between datasets, based on the F1-score.

    Category: benchmark. \n
    Node: critical_difference

    :param metrics_dir_in: The path to the directory, which contains the folder with computed F1-scores.
    :param cd_dir_out: The path to the output directory, which contains the computed results.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_critical_difference_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metrics_dir_in, cd_dir_out)
    return rule
