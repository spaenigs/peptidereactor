import secrets


def _get_header(token):
    return f'''
rule vis_sds_4_Similarity_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(similarity_dir_group_1_in, similarity_dir_group_2_in, html_dir_out):
    return f'''
    input:
         similarity_dir_group_1_in="{similarity_dir_group_1_in}",
         similarity_dir_group_2_in="{similarity_dir_group_2_in}",
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/sds_4_Similarity/Snakefile",
         configfile="nodes/vis/sds_4_Similarity/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(similarity_dir_group_1_in, similarity_dir_group_2_in, html_dir_out,
         benchmark_dir=None):
    """
    Creates a visualization for a single dataset, i.e., a figure composed of
    different heatmaps showing the Phi-correlation- and diversity-based similarity.

    Category: vis \n
    Node: sds_4_Similarity

    :param similarity_dir_group_1_in: The path to the directory pointing to the computed
           similarities for the first group.
    :param similarity_dir_group_2_in: The path to the directory pointing to the computed
           similarities for the second group.
    :param html_dir_out: The path to the output directory containing the Vega-lite specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_sds_4_Similarity_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(similarity_dir_group_1_in, similarity_dir_group_2_in, html_dir_out)
    return rule
