import secrets


def _get_header(token):
    return f'''
rule vis_critical_difference_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(f1_csv_in, crit_diff_dir_in, html_dir_out):
    return f'''
    input:
         f1_csv_in="{f1_csv_in}",
         crit_diff_dir_in="{crit_diff_dir_in}"
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/critical_difference/Snakefile",
         configfile="nodes/vis/critical_difference/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(f1_csv_in, crit_diff_dir_in, html_dir_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_critical_difference_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(f1_csv_in, crit_diff_dir_in, html_dir_out)
    return rule
