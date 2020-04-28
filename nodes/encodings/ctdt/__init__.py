import textwrap


def rule(fasta_in, classes_in, ctdt_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ctdt:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{ctdt_out}"
                params:
                     snakefile="nodes/encodings/ctdt/Snakefile",
                     configfile="nodes/encodings/ctdt/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
