import textwrap


def rule(fasta_in, fasta_out, maps_out):

    rule = textwrap.dedent(
        f'''\
            rule utils_map_sequence_names:
                input:
                     fasta_in="{fasta_in}"
                output:
                     fasta_out="{fasta_out}",
                     maps_out="{maps_out}"
                params:
                     snakefile="nodes/utils/map_sequence_names/Snakefile",
                     configfile="nodes/utils/map_sequence_names/config.yaml"
                run:
                     with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
                         shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")

        ''')

    return rule
