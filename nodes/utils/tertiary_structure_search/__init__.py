import secrets


def _get_header(token):
    return f'''
rule utils_tertiary_structure_search_{token}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{benchmark_out}"'''


def _get_main(fasta_in, classes_in, fasta_sec_out, classes_sec_out, fasta_ter_out, classes_ter_out, pdb_dir, profile_dir):
    return f'''
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
'''


def rule(fasta_in, classes_in, fasta_sec_out, classes_sec_out, fasta_ter_out, classes_ter_out, pdb_dir, profile_dir, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_tertiary_structure_search_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, fasta_sec_out, classes_sec_out, fasta_ter_out, classes_ter_out, pdb_dir, profile_dir)
    return rule

