import textwrap


def rule(fasta_in, classes_in, cgr_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_cgr:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out={cgr_out}
                params:
                     snakefile="nodes/encodings/cgr/Snakefile",
                     configfile="nodes/encodings/cgr/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
