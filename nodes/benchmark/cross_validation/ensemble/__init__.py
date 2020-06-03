import secrets


def _get_header(token):
    return f'''
rule benchmark_cross_validation_ensemble_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(group_1_in, group_2_in, group_1_out, group_2_out):
    return f'''\
    input:
        group_1_in="{group_1_in}",
        group_2_in="{group_2_in}"
    output:
        group_1_out=directory("{group_1_out}"),
        group_2_out=directory("{group_2_out}")
    threads:
         1000      
    params:
        snakefile="nodes/benchmark/cross_validation/ensemble/Snakefile",
        configfile="nodes/benchmark/cross_validation/ensemble/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(group_1_in, group_2_in, group_1_out, group_2_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_cross_validation_ensemble_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(group_1_in, group_2_in, group_1_out, group_2_out)
    return rule
