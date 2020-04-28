import textwrap


def rule(fasta_in, classes_in, gdpc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_gdpc:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{gdpc_out}"
                params:
                     snakefile="nodes/encodings/gdpc/Snakefile",
                     configfile="nodes/encodings/gdpc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
