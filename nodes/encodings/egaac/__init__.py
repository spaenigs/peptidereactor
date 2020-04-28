import textwrap


def rule(fasta_in, classes_in, egaac_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_egaac:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out={egaac_out}
                params:
                     snakefile="nodes/encodings/egaac/Snakefile",
                     configfile="nodes/encodings/egaac/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
