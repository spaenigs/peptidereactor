import secrets


def _get_header(token):
    return f'''
rule encoding_fldpc_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(csv_in, fldpc_out):
    return f'''
    input:
         csv_in="{csv_in}"
    output:
         csv_out={fldpc_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/fldpc/Snakefile",
         configfile="nodes/encodings/fldpc/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(csv_in, fldpc_out, benchmark_dir=None):
    """
    Computes the Five Level Di-Peptide Composition (FLDPC) encoding.

    Category: encodings \n
    Node: fldpc

    :param csv_in: A path pointing to the dataset encoded with the dpc encoding.
    :param fldpc_out: A list of output file paths to store the encoded datasets.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_fldpc_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_in, fldpc_out)
    return rule
