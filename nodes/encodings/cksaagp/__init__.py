import textwrap


def rule(fasta_in, classes_in, length_in, cksaagp_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_cksaagp:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     length_in="{length_in}"
                output:
                     csv_out={cksaagp_out}
                params:
                     snakefile="nodes/encodings/cksaagp/Snakefile",
                     configfile="nodes/encodings/cksaagp/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
