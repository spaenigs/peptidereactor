import textwrap


def rule(csv_in, waac_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_waac:
                input:
                     csv_in="{csv_in}"
                output:
                     csv_out={waac_out}
                params:
                     snakefile="nodes/encodings/waac/Snakefile",
                     configfile="nodes/encodings/waac/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
