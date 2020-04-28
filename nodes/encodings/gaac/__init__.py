import textwrap


def rule(fasta_in, classes_in, gaac_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_gaac:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{gaac_out}"
                params:
                     snakefile="nodes/encodings/gaac/Snakefile",
                     configfile="nodes/encodings/gaac/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
