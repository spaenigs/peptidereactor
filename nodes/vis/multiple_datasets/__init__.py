import secrets


def _get_header(token):
    return f'''
rule vis_multiple_datasets_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(html_files_in, benchmark_csvs_in, html_out):
    return f'''
    input:
         html_files_in={html_files_in},
         benchmark_csv_in={benchmark_csvs_in}
    output:
         html_out="{html_out}"
    threads:
         1000
    params:
         snakefile="nodes/vis/multiple_datasets/Snakefile",
         configfile="nodes/vis/multiple_datasets/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(html_files_in, benchmark_csvs_in, html_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_multiple_datasets_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(html_files_in, benchmark_csvs_in, html_out)
    return rule
