import textwrap


def rule(fasta_in, classes_in, pdb_dir, qsar_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_qsar:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     pdb_dir="{pdb_dir}"
                output:
                     csv_out="{qsar_out}"
                params:
                     snakefile="nodes/encodings/qsar/qsar.smk",
                     configfile="nodes/encodings/qsar/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
