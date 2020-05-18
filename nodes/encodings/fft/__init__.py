import secrets


def _get_header(token):
    return f'''
rule encoding_fft_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(csv_in, fft_out):
    return f'''
    input:
         csv_in={csv_in}
    output:
         csv_out={fft_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/fft/Snakefile",
         configfile="nodes/encodings/fft/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(csv_in, fft_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_fft_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_in, fft_out)
    return rule
