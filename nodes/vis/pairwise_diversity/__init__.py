import secrets


def _get_header(token):
    return f'''
rule vis_pairwise_diversity_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(similarity_dir_group_1_in, similarity_dir_group_2_in,
              ensemble_cv_group_1a_in, ensemble_cv_group_2a_in,
              ensemble_cv_group_1b_in, ensemble_cv_group_2b_in, f1_csv_in, html_dir_out):
    return f'''
    input:
         similarity_dir_group_1_in="{similarity_dir_group_1_in}",
         similarity_dir_group_2_in="{similarity_dir_group_2_in}",
         ensemble_cv_group_1a_in="{ensemble_cv_group_1a_in}", 
         ensemble_cv_group_2a_in="{ensemble_cv_group_2a_in}",
         ensemble_cv_group_1b_in="{ensemble_cv_group_1b_in}", 
         ensemble_cv_group_2b_in="{ensemble_cv_group_2b_in}",
         f1_csv_in="{f1_csv_in}"
    output:
         html_dir_out=directory("{html_dir_out}")
    threads:
         1000
    params:
         snakefile="nodes/vis/pairwise_diversity/Snakefile",
         configfile="nodes/vis/pairwise_diversity/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(similarity_dir_group_1_in, similarity_dir_group_2_in,
         ensemble_cv_group_1a_in, ensemble_cv_group_2a_in,
         ensemble_cv_group_1b_in, ensemble_cv_group_2b_in, f1_csv_in, html_dir_out,
         benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_pairwise_diversity_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(similarity_dir_group_1_in, similarity_dir_group_2_in,
              ensemble_cv_group_1a_in, ensemble_cv_group_2a_in,
              ensemble_cv_group_1b_in, ensemble_cv_group_2b_in, f1_csv_in, html_dir_out)
    return rule
