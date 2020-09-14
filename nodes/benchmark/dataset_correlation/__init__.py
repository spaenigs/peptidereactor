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
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_dataset_correlation_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(group_1_in, group_2_in, metrics_dir_in, dataset_corr_out)
    return rule
