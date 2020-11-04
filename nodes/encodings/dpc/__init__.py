import secrets


def _get_header(token):
    return f'''
rule encoding_dpc_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, dpc_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}"
    output:
         csv_out="{dpc_out}"
    threads:
         1000
    params:
         snakefile="nodes/encodings/dpc/Snakefile",
         configfile="nodes/encodings/dpc/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, dpc_out, benchmark_dir=None):
    """
    Computes the Di-Peptide Composition (DPC) encoding.

    Category: encodings \n
    Node: dpc

    :param fasta_in: The path to the fasta file.
    :param classes_in: The path to the classes file.
    :param dpc_out: The output file path to store the encoded dataset.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_dpc_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, dpc_out)
    return rule
