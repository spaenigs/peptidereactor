import textwrap


def rule(fasta_in, classes_in, profile_dir, disorder_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_disorder:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{disorder_out}"
                params:
                     snakefile="nodes/encodings/disorder/disorder.smk",
                     configfile="nodes/encodings/disorder/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
