import secrets


def _get_header(token):
    return f'''
rule benchmark_similarity_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(group_1_in, group_2_in, corr_dir_out):
    return f'''\
    input:
        group_1_in="{group_1_in}",
        group_2_in="{group_2_in}"
    output:
        corr_dir_out=directory("{corr_dir_out}")
    threads:
         1000      
    params:
        snakefile="nodes/benchmark/similarity/Snakefile",
        configfile="nodes/benchmark/similarity/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(group_1_in, group_2_in, corr_dir_out, benchmark_dir=None):
    """
    Computes the similarity, i.e., the diversity and Phi-correlation of the predictions of
    two groups.

    Category: benchmark. \n
    Node: feature_importance

    :param group_1_in: The path to the directory, which contains the cross-validation results of
           the first group.
    :param group_2_in: The path to the directory, which contains the the cross-validation results of
           the second group.
    :param corr_dir_out: The path to the output directory, to store the computed results.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_similarity_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(group_1_in, group_2_in, corr_dir_out)
    return rule
