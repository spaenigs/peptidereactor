import textwrap


def rule(fasta_in, classes_in, aaindex_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_aaindex:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out={aaindex_out}
                params:
                     snakefile="nodes/encodings/aaindex/Snakefile",
                     configfile="nodes/encodings/aaindex/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
