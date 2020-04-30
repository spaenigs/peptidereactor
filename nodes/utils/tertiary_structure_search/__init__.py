import textwrap


def rule(fasta_in, classes_in, fasta_sec_out, classes_sec_out, fasta_ter_out, classes_ter_out, pdb_dir, profile_dir):

    rule = textwrap.dedent(
        f'''\
            rule util_tertiary_structure_search:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                output:
                     fasta_sec_out="{fasta_sec_out}",
                     classes_sec_out="{classes_sec_out}",
                     profile_dir=directory("{profile_dir}"),
                     fasta_ter_out="{fasta_ter_out}",
                     classes_ter_out="{classes_ter_out}",
                     pdb_dir=directory("{pdb_dir}")
                params:
                     snakefile="nodes/utils/tertiary_structure_search/Snakefile",
                     configfile="nodes/utils/tertiary_structure_search/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
