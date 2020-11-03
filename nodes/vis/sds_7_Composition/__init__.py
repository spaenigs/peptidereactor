import secrets


def _get_header(token):
    return f'''
rule vis_sds_7_Composition_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, html_dir_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}"
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/sds_7_Composition/Snakefile",
         configfile="nodes/vis/sds_7_Composition/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization for a single dataset, i.e., a figure composed of a bar chart
    showing the class distribution, another bar chart depicting the number of sequences
    per length, and a third bar chart visualizing the relative count of amino acids per
    class.

    Category: vis. \n
    Node: sds_7_Composition

    :param fasta_in: A list of fasta file paths.
    :param classes_in: A list of classes file paths.
    :param html_dir_out: The path to the output directory containing the Vega-lite specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_sds_7_Composition_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, html_dir_out)
    return rule
