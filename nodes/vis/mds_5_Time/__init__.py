import secrets


def _get_header(token):
    return f'''
rule vis_mds_5_Time_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fastas_in, benchmark_csvs_in, html_dir_out):
    return f'''
    input:
         fastas_in={fastas_in},
         benchmark_csvs_in={benchmark_csvs_in},    
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/mds_5_Time/Snakefile",
         configfile="nodes/vis/mds_5_Time/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fastas_in, benchmark_csvs_in, html_dir_out, benchmark_dir=None):
    """
    Creates a visualization of multiple datasets, i.e., a figure composed from
    a bar chart (computation time per dataset), a scatter plot (number of sequences vs.
    computation time), and a further scatter plot depicting the average sequence length
    vs. the computation time.

    Category: vis \n
    Node: mds_5_Time

    :param fastas_in: A list of fasta file paths.
    :param benchmark_csvs_in: A list of file paths pointing to benchmarks.
    :param html_dir_out: The path to the output directory containing the Vega-lite specification
           and the accompanied data.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_mds_5_Time_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fastas_in, benchmark_csvs_in, html_dir_out)
    return rule
