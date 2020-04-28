import textwrap


def rule(fasta_in, classes_in, ctriad_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ctriad:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{ctriad_out}"
                params:
                     snakefile="nodes/encodings/ctriad/Snakefile",
                     configfile="nodes/encodings/ctriad/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
