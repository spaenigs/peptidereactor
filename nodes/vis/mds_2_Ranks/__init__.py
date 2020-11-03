import secrets


def _get_header(token):
    return f'''
rule vis_mds_2_Ranks_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(metric_dirs_in, html_dir_out):
    return f'''
    input:
         metric_dirs_in={metric_dirs_in},    
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/mds_2_Ranks/Snakefile",
         configfile="nodes/vis/mds_2_Ranks/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metric_dirs_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization for multiple datasets, i.e., a heatmap sorted by dataset imbalance and
    grouped by encoding type depicting the encoding ranks.

    Category: vis \n
    Node: mds_2_Ranks

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
        benchmark_out = f"{benchmark_dir}vis_mds_2_Ranks_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metric_dirs_in, html_dir_out)
    return rule
