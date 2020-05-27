import secrets


def _get_header(token):
    return f'''
rule utils_check_dataset_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, report_out):
    return f'''
    input:
         fasta_in="{fasta_in}"
    output:
         report_out="{report_out}"
    threads:
         1000
    priority:
         1000 
    params:
         snakefile="nodes/utils/check_dataset/Snakefile",
         configfile="nodes/utils/check_dataset/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, report_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_check_dataset_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, report_out)
    return rule
