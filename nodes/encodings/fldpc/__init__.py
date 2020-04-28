import textwrap


def rule(csv_in, fldpc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_fldpc:
                input:
                     csv_in="{csv_in}"
                output:
                     csv_out={fldpc_out}
                params:
                     snakefile="nodes/encodings/fldpc/Snakefile",
                     configfile="nodes/encodings/fldpc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
