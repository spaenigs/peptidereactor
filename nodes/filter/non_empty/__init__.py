import textwrap


def rule(csv_in, csv_out):

    rule = textwrap.dedent(
        f'''\
            rule filter_non_empty:
                input:
                     csv_dir_in={csv_in}
                output:
                     csv_out=directory("{csv_out}")
                params:
                     snakefile="nodes/filter/non_empty/Snakefile",
                     configfile="nodes/filter/non_empty/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

    ''')

    return rule