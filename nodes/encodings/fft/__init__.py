import textwrap


def rule(csv_in, fft_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_fft:
                input:
                     csv_in={csv_in}
                output:
                     csv_out={fft_out}
                params:
                     snakefile="nodes/encodings/fft/Snakefile",
                     configfile="nodes/encodings/fft/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
