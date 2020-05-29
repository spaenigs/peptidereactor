import secrets


def _get_header(token):
    return f'''
rule benchmark_feature_importance_{token}:
'''


def _get_benchmark(benchmark_out):
    return f'''\
    benchmark:
        "{benchmark_out}"
'''


def _get_main(feat_imp_in, feat_imp_out):
    return f'''\
    input:
        feat_imp_in="{feat_imp_in}"
    output:
        feat_imp_out="{feat_imp_out}"
    threads:
         1000
    params:
        snakefile="nodes/benchmark/feature_importance/Snakefile",
        configfile="nodes/benchmark/feature_importance/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(feat_imp_in, feat_imp_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}benchmark_feature_importance_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(feat_imp_in, feat_imp_out)
    return rule
