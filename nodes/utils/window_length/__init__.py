import secrets


def _get_header(token):
    return f'''
rule utils_window_length_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, length_out):
    return f'''
    input:
         fasta_in="{fasta_in}"
    output:
         length_out={length_out}
    threads:
         1000
    params:
         snakefile="nodes/utils/window_length/Snakefile",
         configfile="nodes/utils/window_length/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, length_out, benchmark_dir=None):
    """
    Computes the maximum allowed values for parameter-based encodings.

    Category: utils. \n
    Node: window_length

    :param fasta_in: The path to the fasta file.
    :param length_out: The path to the output file for a specific encoding.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_window_length_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, length_out)
    return rule
