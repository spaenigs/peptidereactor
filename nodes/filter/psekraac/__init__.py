import secrets


def _get_header(token):
    return f'''
rule filter_psekraac_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(csv_in, csv_out):
    return f'''
    input:
         csv_in={csv_in}
    output:
         csv_out=directory("{csv_out}")
    params:
         snakefile="nodes/filter/psekraac/Snakefile",
         configfile="nodes/filter/psekraac/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(csv_in, csv_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}filter_psekraac_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(csv_in, csv_out)
    return rule
