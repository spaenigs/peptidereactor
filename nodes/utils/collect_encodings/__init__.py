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
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_collect_encodings_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_seq_in, csv_str_in, csv_seq_out, csv_str_out)
    return rule
