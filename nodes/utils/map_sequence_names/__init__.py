import secrets


def _get_header(token):
    return f'''
rule utils_map_sequence_names_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, fasta_out, maps_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}"
    output:
         fasta_out="{fasta_out}",
         maps_out="{maps_out}"
    threads:
         1000
    params:
         snakefile="nodes/utils/map_sequence_names/Snakefile",
         configfile="nodes/utils/map_sequence_names/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, fasta_out, classes_in, maps_out, benchmark_dir=None):
    """
    Maps the names of the input sequences to a general naming structure. In addition,
    this node ensures a valid input format, e.g., valid amino acid sequences.

    Category: utils. \n
    Node: map_sequence_names

    :param fasta_in: The path to the fasta file.
    :param classes_in: The path to the classes file.
    :param fasta_out: The output path to the sanitized fasta file.
    :param maps_out: The path to the output file containing the mappings from the original names
           to the new names.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_map_sequence_names_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, fasta_out, maps_out)
    return rule

