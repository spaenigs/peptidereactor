import secrets


def _get_header(token):
    return f'''
rule encoding_psekraac_type16_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, psekraac_type16_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}"
    output:
         csv_out={psekraac_type16_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/psekraac_type16/Snakefile",
         configfile="nodes/encodings/psekraac_type16/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, psekraac_type16_out, benchmark_dir=None):
    """
    Computes the pseudo K-tuple reduced amino acids composition (PseKRAAC)
    type 16 encoding.

    Category: encodings \n
    Node: psekraac_type16

    :param fasta_in: The path to the fasta file.
    :param classes_in: The path to the classes file.
    :param psekraac_type16_out: A list of output file paths to store the encoded datasets.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_psekraac_type16_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, psekraac_type16_out)
    return rule
