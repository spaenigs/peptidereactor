import textwrap


def rule(fasta_in, classes_in, fasta_out):

    rule = textwrap.dedent(
        f'''\
            rule util_multiple_sequence_alignment:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     fasta_out="{fasta_out}"
                params:
                     snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
                     configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

        ''')

    return rule