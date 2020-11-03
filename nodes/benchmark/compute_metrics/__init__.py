import secrets


def _get_header(token):
    return f'''
rule benchmark_compute_metrics_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(csv_dir_in, metrics_dir_out):
    return f'''\
    input:
        csv_dir_in="{csv_dir_in}"
    output:
        metrics_dir_out=directory("{metrics_dir_out}")
    threads:
         1000
    params:
        snakefile="nodes/benchmark/compute_metrics/Snakefile",
        configfile="nodes/benchmark/compute_metrics/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(csv_dir_in, metrics_dir_out, benchmark_dir=None):
    """
    Computes the F1-score, MCC, Precision, Recall, Sensitivity, and Specificity
    for given cross-validation results.

    Category: benchmark. \n
    Node: compute_metrics

    :param csv_dir_in: The path to the directory containing the matrices of size splits x fold_size.
    :param metrics_dir_out: The path to the output directory, which contains the computed metrics.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_compute_metrics_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_dir_in, metrics_dir_out)
    return rule
