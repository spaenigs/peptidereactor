import textwrap


def rule(ngram_type,
        fasta_in=None, classes_in=None, csv_in=None, length_in=None,
        ngram_out=None, ngram_lsv_out=None, ngram_sv_out=None):

    if ngram_type in ["a2", "a3"]:
        rule = textwrap.dedent(
            f'''\
                rule encoding_ngram_{ngram_type}:
                    input:
                         csv_in="{csv_in}",
                         length_in="{length_in}"
                    output:
                         csv_out={ngram_out},
                         lsv_out={ngram_lsv_out},
                         sv_out={ngram_sv_out}
                    params:
                         snakefile="nodes/encodings/ngram/Snakefile",
                         configfile="nodes/encodings/ngram/config.yaml"
                    run:
                         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

                ''')

    elif ngram_type in ["e2", "e3", "s2", "s3"]:
        rule = textwrap.dedent(
            f'''\
                rule encoding_ngram_{ngram_type}:
                    input:
                         fasta_in="{fasta_in}",
                         classes_in="{classes_in}",
                         length_in="{length_in}"
                    output:
                         csv_out={ngram_out},
                         lsv_out={ngram_lsv_out},
                         sv_out={ngram_sv_out}
                    params:
                         snakefile="nodes/encodings/ngram/Snakefile",
                         configfile="nodes/encodings/ngram/config.yaml"
                    run:
                         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                             shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

                ''')

    else:
        raise ValueError(f"Unknown type for ngram encoding: {ngram_type}")

    return rule
