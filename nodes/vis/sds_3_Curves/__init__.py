import secrets


def _get_header(token):
    return f'''
rule vis_sds_3_Curves_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(metrics_dir_in, html_dir_out):
    return f'''
    input:
         metrics_dir_in="{metrics_dir_in}",
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/sds_3_Curves/Snakefile",
         configfile="nodes/vis/sds_3_Curves/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metrics_dir_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization for a single dataset, i.e., a figure composed of
    the ROC- as well as the Precision-Recall-curve for the top-6-encodings and
    the top-3-sequence- and top-3-structure-based encodings.

    Category: vis \n
    Node: sds_3_Curves

    :param metric_dirs_in: A list of directory paths pointing to the metric directory.
    :param html_dir_out: The path to the output directory containing the Vega-lite specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_sds_3_Curves_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metrics_dir_in, html_dir_out)
    return rule
