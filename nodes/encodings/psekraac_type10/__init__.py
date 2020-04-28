import textwrap


def rule(fasta_in, classes_in, psekraac_type10_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_psekraac_type10:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out={psekraac_type10_out}
                params:
                     snakefile="nodes/encodings/psekraac_type10/Snakefile",
                     configfile="nodes/encodings/psekraac_type10/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
