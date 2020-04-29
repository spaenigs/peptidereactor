import textwrap


def rule(fasta_in, classes_in, profile_dir, ta_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ta:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{ta_out}"
                params:
                     snakefile="nodes/encodings/ta/ta.smk",
                     configfile="nodes/encodings/ta/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
