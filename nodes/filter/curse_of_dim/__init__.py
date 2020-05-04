import textwrap


def rule(csv_in, csv_out):

    rule = textwrap.dedent(
        f'''\
            rule filter_curse_of_dim:
                input:
                     csv_dir_in={csv_in}
                output:
                     csv_out=directory("{csv_out}")
                params:
                     snakefile="nodes/filter/curse_of_dim/Snakefile",
                     configfile="nodes/filter/curse_of_dim/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

    ''')

    return rule