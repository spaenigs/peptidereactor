import secrets


def _get_header(token):
    return f'''
rule encoding_gdpc_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, gdpc_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}"
    output:
         csv_out="{gdpc_out}"
    threads:
         1000
    params:
         snakefile="nodes/encodings/gdpc/Snakefile",
         configfile="nodes/encodings/gdpc/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, gdpc_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_gdpc_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, gdpc_out)
    return rule
