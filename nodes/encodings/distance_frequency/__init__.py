import textwrap


def rule(fasta_in, classes_in, distance_frequency_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_distance_frequency:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out={distance_frequency_out}
                params:
                     snakefile="nodes/encodings/distance_frequency/Snakefile",
                     configfile="nodes/encodings/distance_frequency/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
