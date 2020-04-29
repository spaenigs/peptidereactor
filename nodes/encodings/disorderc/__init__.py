import textwrap


def rule(fasta_in, classes_in, profile_dir, disorderc_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_disorderc:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{disorderc_out}"
                params:
                     snakefile="nodes/encodings/disorderc/disorderc.smk",
                     configfile="nodes/encodings/disorderc/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule