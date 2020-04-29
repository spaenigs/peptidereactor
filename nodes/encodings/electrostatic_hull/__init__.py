import textwrap


def rule(fasta_in, classes_in, pdb_dir, electrostatic_hull_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_electrostatic_hull:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     pdb_dir="{pdb_dir}"
                output:
                     csv_out={electrostatic_hull_out}
                params:
                     snakefile="nodes/encodings/electrostatic_hull/electrostatic_hull.smk",
                     configfile="nodes/encodings/electrostatic_hull/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
