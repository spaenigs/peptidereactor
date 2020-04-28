import textwrap


def rule(fasta_in, classes_in, aac_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_aac:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{aac_out}"
                params:
                     snakefile="nodes/encodings/aac/Snakefile",
                     configfile="nodes/encodings/aac/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
