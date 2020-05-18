import secrets


def _get_header(token, ngram_type):
    return f'''
rule encoding_ngram_{ngram_type}_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(ngram_type,
              fasta_in=None, classes_in=None, csv_in=None, length_in=None,
              ngram_out=None, ngram_lsv_out=None, ngram_sv_out=None):

    if ngram_type in ["a2", "a3"]:
        return f'''
    input:
         csv_in="{csv_in}",
         length_in="{length_in}"
    output:
         csv_out={ngram_out},
         lsv_out={ngram_lsv_out},
         sv_out={ngram_sv_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"'''

    elif ngram_type in ["e2", "e3", "s2", "s3"]:
        return f'''
    input:
         fasta_in="{fasta_in}",
         classes_in="{classes_in}",
         length_in="{length_in}"
    output:
         csv_out={ngram_out},
         lsv_out={ngram_lsv_out},
         sv_out={ngram_sv_out}
    threads:
         1000
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"'''

    else:
        raise ValueError(f"Unknown type for ngram encoding: {ngram_type}")


def _get_footer():
    return f'''
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(ngram_type,
         fasta_in=None, classes_in=None, csv_in=None, length_in=None,
         ngram_out=None, ngram_lsv_out=None, ngram_sv_out=None,
         benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token, ngram_type)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_ngram_{ngram_type}_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(ngram_type,
                      fasta_in=fasta_in, classes_in=classes_in, csv_in=csv_in, length_in=length_in,
                      ngram_out=ngram_out, ngram_lsv_out=ngram_lsv_out, ngram_sv_out=ngram_sv_out)
    rule += _get_footer()
    return rule
