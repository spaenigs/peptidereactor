import textwrap


def rule(fasta_in, classes_in, ctdd_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ctdd:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{ctdd_out}"
                params:
                     snakefile="nodes/encodings/ctdd/Snakefile",
                     configfile="nodes/encodings/ctdd/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
