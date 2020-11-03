import secrets


def _get_header(token):
    return f'''
rule utils_collect_encodings_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(csv_seq_in, csv_str_in, csv_seq_out, csv_str_out):
    return f'''
    input:
         csv_seq_in={csv_seq_in},
         csv_str_in={csv_str_in}
    output:
         csv_seq_out=directory("{csv_seq_out}"),
         csv_str_out=directory("{csv_str_out}")
    threads:
         1000
    params:
         snakefile="nodes/utils/collect_encodings/Snakefile",
         configfile="nodes/utils/collect_encodings/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(csv_seq_in, csv_str_in, csv_seq_out, csv_str_out, benchmark_dir=None):
    """
    Collects the encoded datasets into the specified directory.

    Category: utils. \n
    Node: collect_encodings

    :param csv_seq_in: A list of output file paths pointing sequence-based datasets.
    :param csv_str_in: A list of output file paths pointing structure-based datasets.
    :param csv_seq_out: The path to the output file.
    :param csv_str_out: The path to the output file.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).

    :return: A string object representing a Snakemake rule.
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_collect_encodings_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_seq_in, csv_str_in, csv_seq_out, csv_str_out)
    return rule
