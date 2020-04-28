import textwrap


def rule(fasta_in, classes_in, blomap_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_blomap:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{blomap_out}"
                params:
                     snakefile="nodes/encodings/blomap/Snakefile",
                     configfile="nodes/encodings/blomap/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
