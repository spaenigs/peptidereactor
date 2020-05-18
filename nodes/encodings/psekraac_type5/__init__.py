import secrets


def _get_header(token):
    return f'''
rule encoding_psekraac_type5_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, psekraac_type5_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}"
    output:
         csv_out={psekraac_type5_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/psekraac_type5/Snakefile",
         configfile="nodes/encodings/psekraac_type5/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, psekraac_type5_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_psekraac_type5_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, psekraac_type5_out)
    return rule
