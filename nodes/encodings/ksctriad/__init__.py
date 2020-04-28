import textwrap


def rule(fasta_in, classes_in, length_in, ksctriad_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ksctriad:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     length_in="{length_in}"
                output:
                     csv_out={ksctriad_out}
                params:
                     snakefile="nodes/encodings/ksctriad/Snakefile",
                     configfile="nodes/encodings/ksctriad/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
