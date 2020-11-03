import secrets


def _get_header(token):
    return f'''
rule benchmark_dataset_correlation_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(group_1_in, group_2_in, metrics_dir_in, dataset_corr_out):
    return f'''\
    input:
        group_1_in="{group_1_in}",
        group_2_in="{group_2_in}",
        metrics_dir_in="{metrics_dir_in}"
    output:
        dataset_corr_out="{dataset_corr_out}"
    threads:
         1000      
    params:
        snakefile="nodes/benchmark/dataset_correlation/Snakefile",
        configfile="nodes/benchmark/dataset_correlation/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(group_1_in, group_2_in, metrics_dir_in, dataset_corr_out, benchmark_dir=None):
    """
    Computes the dataset correlation of the best datasets of two groups, determined by the F1-score.

    Category: benchmark. \n
    Node: dataset_correlation

    :param group_1_in: The path to the directory, which contains the encoded datasets of the first group.
    :param group_2_in: The path to the directory, which contains the encoded datasets of the second group.
    :param metrics_dir_in: The path to the directory, which contains the folder with computed F1-scores.
    :param dataset_corr_out: The path to the output directory, to store the computed results.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_dataset_correlation_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(group_1_in, group_2_in, metrics_dir_in, dataset_corr_out)
    return rule
