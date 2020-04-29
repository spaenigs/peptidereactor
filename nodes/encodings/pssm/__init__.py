import textwrap


def rule(fasta_in, classes_in, profile_dir, pssm_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_pssm:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     profile="{profile_dir}"
                output:
                     csv_out="{pssm_out}"
                params:
                     snakefile="nodes/encodings/pssm/Snakefile",
                     configfile="nodes/encodings/pssm/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
