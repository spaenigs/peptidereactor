import secrets


def _get_header(token):
    return f'''
rule vis_sds_6_Difference_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(metrics_dir_in, crit_diff_dir_in, html_dir_out):
    return f'''
    input:
         metrics_dir_in="{metrics_dir_in}",
         crit_diff_dir_in="{crit_diff_dir_in}"
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/sds_6_Difference/Snakefile",
         configfile="nodes/vis/sds_6_Difference/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metrics_dir_in, crit_diff_dir_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization for a single dataset, i.e., a figure composed of a
    heatmap showing the critical differences for all vs. all encodings, a bar chart
    depicting the number of critical different datasets per encoding group, and the
    binned number of critical different encodings.

    Category: vis \n
    Node: sds_6_Difference

    :param metrics_dir_in: The path to the metrics directory.
    :param crit_diff_dir_in: The path to the results of the critical differences.
    :param html_dir_out: The path to the output directory containing the Vega-lite specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_sds_6_Difference_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metrics_dir_in, crit_diff_dir_in, html_dir_out)
    return rule
