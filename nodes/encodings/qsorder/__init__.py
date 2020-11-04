import secrets


def _get_header(token):
    return f'''
rule encoding_qsorder_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, length_in, qsorder_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}",
         length_in="{length_in}"
    output:
         csv_out={qsorder_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/qsorder/Snakefile",
         configfile="nodes/encodings/qsorder/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, length_in, qsorder_out, benchmark_dir=None):
    """
    Computes the Quasi-sequence-order (QSOrder) encoding.

    Category: encodings \n
    Node: qsorder

    :param fasta_in: The path to the fasta file.
    :param classes_in: The path to the classes file.
    :param length_in: The path to the file, containing the allowed parameter space.
    :param qsorder_out: A list of output file paths to store the encoded datasets.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_qsorder_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, length_in, qsorder_out)
    return rule
