import secrets


def _get_header(token):
    return f'''
rule vis_md_ranks_hm_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(metric_dirs_in, html_dir_out):
    return f'''
    input:
         metric_dirs_in={metric_dirs_in},    
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/md_ranks_hm/Snakefile",
         configfile="nodes/vis/md_ranks_hm/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(metric_dirs_in, html_dir_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_md_ranks_hm_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(metric_dirs_in, html_dir_out)
    return rule
