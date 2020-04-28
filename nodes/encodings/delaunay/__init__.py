import textwrap


def rule(fasta_in, classes_in, pdb_dir, delaunay_out):

    rule = textwrap.dedent(
        f'''\
            rule encoding_delaunay:
                input:
                     fasta_in="{fasta_in}",
                     classes_in="{classes_in}",
                     pdb_dir="{pdb_dir}"
                output:
                     csv_out={delaunay_out}
                params:
                     snakefile="nodes/encodings/delaunay/delaunay.smk",
                     configfile="nodes/encodings/delaunay/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

            ''')

    return rule
