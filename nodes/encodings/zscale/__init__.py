import textwrap


def rule(fasta_in, classes_in, zscale_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_zscale:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{zscale_out}"
                params:
                     snakefile="nodes/encodings/zscale/Snakefile",
                     configfile="nodes/encodings/zscale/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
