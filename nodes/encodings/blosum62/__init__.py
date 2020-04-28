import textwrap


def rule(fasta_in, classes_in, blosum62_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_blosum62:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}"
                output:
                     csv_out="{blosum62_out}"
                params:
                     snakefile="nodes/encodings/blosum62/Snakefile",
                     configfile="nodes/encodings/blosum62/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
