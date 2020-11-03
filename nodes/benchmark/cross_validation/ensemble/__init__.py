import secrets


def _get_header(token):
    return f'''
rule benchmark_cross_validation_ensemble_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(group_1_in, group_2_in, group_1_out, group_2_out):
    return f'''\
    input:
        group_1_in="{group_1_in}",
        group_2_in="{group_2_in}"
    output:
        group_1_out=directory("{group_1_out}"),
        group_2_out=directory("{group_2_out}")
    threads:
         1000      
    params:
        snakefile="nodes/benchmark/cross_validation/ensemble/Snakefile",
        configfile="nodes/benchmark/cross_validation/ensemble/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(group_1_in, group_2_in, group_1_out, group_2_out, benchmark_dir=None):
    """
    Executes cross validation between all datasets of two groups. For later usage as part of an ensemble
    classifier, this node ensures that all datasets contain the same elements.

    Category: benchmark. \n
    Node: ensemble

    :param group_1_in: The path to the directory, which contains the encoded datasets of the first group.
    :param group_2_in: The path to the directory, which contains the encoded datasets of the second group.
    :param group_1_out: The path to the output directory, to store the computed results for the first group.
    :param group_2_out: The path to the output directory, to store the computed results for the second group.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_cross_validation_ensemble_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(group_1_in, group_2_in, group_1_out, group_2_out)
    return rule
