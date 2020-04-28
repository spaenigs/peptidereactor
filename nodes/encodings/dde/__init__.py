import textwrap


def rule(fasta_in, classes_in, dde_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_dde:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{dde_out}"
                params:
                     snakefile="nodes/encodings/dde/Snakefile",
                     configfile="nodes/encodings/dde/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
