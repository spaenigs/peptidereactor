import secrets


def _get_header(token):
    return f'''
rule vis_single_dataset_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, encoding_benchmark_dir_in, benchmark_csv_in, html_out):
    return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}",
         metrics_dir_in="{encoding_benchmark_dir_in + "metrics/"}",
         dataset_correlation_in="{encoding_benchmark_dir_in + "dataset_correlation.csv"}",
         similarity_dir_group_1_in="{encoding_benchmark_dir_in + "similarity/seq_vs_str/"}",
         similarity_dir_group_2_in="{encoding_benchmark_dir_in + "similarity/all_vs_all/"}",
         ensemble_cv_group_1a_in="{encoding_benchmark_dir_in + "ensemble/seq_vs_str/sequence_based/"}", 
         ensemble_cv_group_2a_in="{encoding_benchmark_dir_in + "ensemble/seq_vs_str/structure_based/"}",
         ensemble_cv_group_1b_in="{encoding_benchmark_dir_in + "ensemble/all_vs_all/group_1/"}", 
         ensemble_cv_group_2b_in="{encoding_benchmark_dir_in + "ensemble/all_vs_all/group_2/"}",
         crit_diff_dir_in="{encoding_benchmark_dir_in + "friedman/"}",
         benchmark_csv_in="{benchmark_csv_in}"
    output:
         html_out="{html_out}"
    threads:
         1000
    params:
         snakefile="nodes/vis/single_dataset/Snakefile",
         configfile="nodes/vis/single_dataset/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, encoding_benchmark_dir_in, benchmark_csv_in, html_out,
         benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}vis_single_dataset_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, encoding_benchmark_dir_in, benchmark_csv_in, html_out)
    return rule
