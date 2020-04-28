import textwrap


def rule(fasta_in, classes_in, length_in, qsorder_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_qsorder:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     length_in="{length_in}"
                output:
                     csv_out={qsorder_out}
                params:
                     snakefile="nodes/encodings/qsorder/Snakefile",
                     configfile="nodes/encodings/qsorder/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
