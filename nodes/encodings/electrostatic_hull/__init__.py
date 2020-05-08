import secrets


def _get_header(token):
    return f'''
rule encoding_electrostatic_hull_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, pdb_dir, electrostatic_hull_out):
    return f'''
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
'''


def rule(fasta_in, classes_in, pdb_dir, electrostatic_hull_out, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}encoding_electrostatic_hull_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, pdb_dir, electrostatic_hull_out)
    return rule
