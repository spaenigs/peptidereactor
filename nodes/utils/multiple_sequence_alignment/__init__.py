import textwrap


def rule(fastas_in, fastas_out):

    rule = textwrap.dedent(
        f'''\
            rule util_multiple_sequence_alignment:
                input:
                     fastas_in={fastas_in}
                output:
                     fastas_out={fastas_out}
                params:
                     snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
                     configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

        ''')

    return rule