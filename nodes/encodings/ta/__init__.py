import secrets


def _get_header(token):
    return f'''
rule encoding_ta_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, profile_dir, ta_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}",
         profile="{profile_dir}"
    output:
         csv_out="{ta_out}"
    threads:
         1000
    params:
         snakefile="nodes/encodings/ta/ta.smk",
         configfile="nodes/encodings/ta/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, profile_dir, ta_out, benchmark_dir=None):
    """
    Computes the Torsion angle (TA) encoding. This encoding is based on
    the secondary structure.

    Category: encodings \n
    Node: ta

    :param fasta_in: The path to the fasta file.
    :param classes_in: The path to the classes file.
    :param profile_dir: The path to the directory, containing the pre-computed profiles.
    :param ta_out: The output file path to store the encoded dataset.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_ta_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, profile_dir, ta_out)
    return rule
