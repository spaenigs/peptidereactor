import textwrap


def rule(fasta_in, classes_in, tpc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_tpc:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{tpc_out}"
                params:
                     snakefile="nodes/encodings/tpc/Snakefile",
                     configfile="nodes/encodings/tpc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
