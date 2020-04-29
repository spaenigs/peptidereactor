import textwrap


def rule(fasta_in, classes_in, profile_dir, sseb_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_sseb:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{sseb_out}"
                params:
                     snakefile="nodes/encodings/sseb/sseb.smk",
                     configfile="nodes/encodings/sseb/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
