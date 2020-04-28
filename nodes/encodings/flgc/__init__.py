import textwrap


def rule(csv_in, flgc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_flgc:
                input:
                     csv_in="{csv_in}"
                output:
                     csv_out={flgc_out}
                params:
                     snakefile="nodes/encodings/flgc/Snakefile",
                     configfile="nodes/encodings/flgc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
