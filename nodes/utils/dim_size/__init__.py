import textwrap


def rule(fasta_in, length_out):

    rule = textwrap.dedent(
        f'''\
            rule util_dim_size:
                input:
                     fasta_in="{fasta_in}"
                output:
                     length_out={length_out}
                params:
                     snakefile="nodes/utils/dim_size/Snakefile",
                     configfile="nodes/utils/dim_size/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

    ''')

    return rule