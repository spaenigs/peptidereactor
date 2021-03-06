import secrets


def _get_header(token):
    return f'''
rule vis_sds_8_Correlation_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(metrics_dir_in, dataset_correlation_in, html_dir_out):
    return f'''
    input:
         metrics_dir_in="{metrics_dir_in}",
         dataset_correlation_in="{dataset_correlation_in}"
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/sds_8_Correlation/Snakefile",
         configfile="nodes/vis/sds_8_Correlation/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metrics_dir_in, dataset_correlation_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization for a single dataset, i.e., a clustered dendrogram of
    correlated datasets.

    Category: vis \n
    Node: sds_8_Correlation

    :param metrics_dir_in: A list of directory paths pointing to the metric directory.
    :param dataset_correlation_in: The path to the results of the dataset correlation.
    :param html_dir_out: The path to the output directory containing the Vega specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_sds_8_Correlation_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metrics_dir_in, dataset_correlation_in, html_dir_out)
    return rule
