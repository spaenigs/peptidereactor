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


def _get_main(csv_seq_in, csv_str_in, csv_dir_out):
    return f'''\
    input:
        csv_seq_in="{csv_seq_in}",
        csv_str_in="{csv_str_in}"
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


def rule(csv_seq_in, csv_str_in, csv_dir_out, benchmark_dir=None):
    """
    Executes cross validation between all datasets of two groups.

    Category: benchmark. \n
    Node: single

    :param csv_seq_in: The path to the directory, which contains the encoded datasets of the first group.
    :param csv_str_in: The path to the directory, which contains the encoded datasets of the second group.
    :param csv_dir_out: The path to the output directory, to store the computed results.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_cross_validation_single_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_seq_in, csv_str_in, csv_dir_out)
    return rule
