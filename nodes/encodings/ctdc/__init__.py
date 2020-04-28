import textwrap


def rule(fasta_in, classes_in, ctdc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ctdc:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{ctdc_out}"
                params:
                     snakefile="nodes/encodings/ctdc/Snakefile",
                     configfile="nodes/encodings/ctdc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
