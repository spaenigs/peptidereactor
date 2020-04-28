import textwrap


def rule(fasta_in, classes_in, dpc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_dpc:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{dpc_out}"
                params:
                     snakefile="nodes/encodings/dpc/Snakefile",
                     configfile="nodes/encodings/dpc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
