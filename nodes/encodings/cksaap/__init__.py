import textwrap


def rule(fasta_in, classes_in, length_in, cksaap_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_cksaap:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     length_in="{length_in}"
                output:
                     csv_out={cksaap_out}
                params:
                     snakefile="nodes/encodings/cksaap/Snakefile",
                     configfile="nodes/encodings/cksaap/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
