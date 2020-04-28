import textwrap


def rule(fasta_in, classes_in, length_in, socnumber_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_socnumber:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     length_in="{length_in}"
                output:
                     csv_out={socnumber_out}
                params:
                     snakefile="nodes/encodings/socnumber/Snakefile",
                     configfile="nodes/encodings/socnumber/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
