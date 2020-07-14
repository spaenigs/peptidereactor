import secrets


def _get_header(token):
    return f'''
rule utils_collect_benchmark_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(final_dirs_in, final_files_in, csv_out):
    return f'''
    input:
         final_dirs_in={final_dirs_in},
         final_files_in={final_files_in}
    output:
         csv_out="{csv_out}"
    threads:
         1000
    params:
         snakefile="nodes/utils/collect_benchmark/Snakefile",
         configfile="nodes/utils/collect_benchmark/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(final_dirs_in, final_files_in, csv_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_collect_benchmark_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(final_dirs_in, final_files_in, csv_out)
    return rule
