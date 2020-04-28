import textwrap


def rule(fasta_in, classes_in, binary_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_binary:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{binary_out}"
                params:
                     snakefile="nodes/encodings/binary/Snakefile",
                     configfile="nodes/encodings/binary/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
