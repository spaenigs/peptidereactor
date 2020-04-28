import textwrap


def rule(fasta_in, classes_in, length_in, paac_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_paac:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     length_in="{length_in}"
                output:
                     csv_out={paac_out}
                params:
                     snakefile="nodes/encodings/paac/Snakefile",
                     configfile="nodes/encodings/paac/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
