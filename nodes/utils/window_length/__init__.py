import textwrap


def rule(fasta_in, length_out):

    rule = textwrap.dedent(
        f'''\
            rule util_window_length:
                input:
                     fasta_in="{fasta_in}"
                output:
                     length_out={length_out}
                params:
                     snakefile="nodes/utils/window_length/Snakefile",
                     configfile="nodes/utils/window_length/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

    ''')

    return rule