import textwrap


def rule(fasta_in, classes_in, profile_dir, ssec_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_ssec:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{ssec_out}"
                params:
                     snakefile="nodes/encodings/ssec/ssec.smk",
                     configfile="nodes/encodings/ssec/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
