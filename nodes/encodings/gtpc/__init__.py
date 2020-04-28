import textwrap


def rule(fasta_in, classes_in, gtpc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_gtpc:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{gtpc_out}"
                params:
                     snakefile="nodes/encodings/gtpc/Snakefile",
                     configfile="nodes/encodings/gtpc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
