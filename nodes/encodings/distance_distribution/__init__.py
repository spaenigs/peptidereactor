import textwrap


def rule(fasta_in, classes_in, pdb_dir, distance_distribution_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_distance_distribution:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     pdb_dir="{pdb_dir}"
                output:
                     csv_out="{distance_distribution_out}"
                params:
                     snakefile="nodes/encodings/distance_distribution/distance_distribution.smk",
                     configfile="nodes/encodings/distance_distribution/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
