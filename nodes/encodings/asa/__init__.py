import textwrap


def rule(fasta_in, classes_in, profile_dir, asa_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_asa:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{asa_out}"
                params:
                     snakefile="nodes/encodings/asa/asa.smk",
                     configfile="nodes/encodings/asa/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
