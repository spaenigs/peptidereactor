import secrets


def _get_header(token):
    return f'''
rule utils_multiple_sequence_alignment_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fastas_in, fastas_out):
    return f'''
    input:
         fastas_in={fastas_in}
    output:
         fastas_out={fastas_out}
    threads:
         1000
    params:
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fastas_in, fastas_out, benchmark_dir=None):
    """
    Conducts the multiple sequence alignment.

    Category: utils. \n
    Node: multiple_sequence_alignment

    :param fastas_in: A list of file paths pointing to fasta files.
    :param fastas_out: A list of output file paths.
    :param maps_out: The path to the output file containing the mappings from the original names
           to the new names.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_multiple_sequence_alignment_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fastas_in, fastas_out)
    return rule

