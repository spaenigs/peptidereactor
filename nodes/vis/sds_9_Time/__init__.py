import secrets


def _get_header(token):
    return f'''
rule vis_sds_9_Time_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(metrics_dir_in, benchmark_csv_in, html_dir_out):
    return f'''
    input:
         metrics_dir_in="{metrics_dir_in}",
         benchmark_csv_in="{benchmark_csv_in}"
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/sds_9_Time/Snakefile",
         configfile="nodes/vis/sds_9_Time/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metrics_dir_in, benchmark_csv_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization for a single dataset, i.e., a figure composed of a
    scatter plot, showing the median performance of the encoding groups and
    a further scatter plot, depicting the required computation time for a specific
    node.

    Category: vis \n
    Node: sds_9_Time

    :param metrics_dir_in: The path to the metrics directory.
    :param benchmark_csv_in: The path to the benchmark file.
    :param html_dir_out: The path to the output directory containing the Vega-lite specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_sds_9_Time_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metrics_dir_in, benchmark_csv_in, html_dir_out)
    return rule
