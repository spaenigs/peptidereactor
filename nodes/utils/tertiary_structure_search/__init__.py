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
    threads:
         1000
    params:
         snakefile="nodes/utils/tertiary_structure_search/Snakefile",
         configfile="nodes/utils/tertiary_structure_search/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f"""{{e.snakemake}} -s {{params.snakefile}} --configfile {{params.configfile}}""")
'''


def rule(fasta_in, classes_in, fasta_sec_out, classes_sec_out,
         fasta_ter_out, classes_ter_out, pdb_dir, profile_dir,
         benchmark_dir=None):
    """
    Approximates the tertiary structure of the input sequences. Note that, not for
    all sequences a structure can be found.

    Category: utils. \n
    Node: tertiary_structure_search

    :param fasta_in: The path to the fasta file.
    :param classes_in: The path to the classes file.
    :param fasta_sec_out: The output path to the fasta file containing the sequences
           for which a secondary structure could have been found.
    :param classes_sec_out: The output path to the classes file containing the classes of
           sequences for which a secondary structure could have been found.
    :param fasta_ter_out: The output path to the fasta file containing the sequences
           for which a tertiary structure could have been found.
    :param classes_ter_out: The output path to the classes file containing the classes of
           sequences for which a secondary structure could have been found.
    :param pdb_dir: The path to the output directory to store the structures.
    :param profile_dir: The path to the output directory to store the profiles.
    :param benchmark_dir: The path to the directory to store the benchmark results. If None,
           benchmark will be not executed (default).
    """
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{benchmark_dir}utils_tertiary_structure_search_{token}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main(fasta_in, classes_in, fasta_sec_out, classes_sec_out, fasta_ter_out, classes_ter_out, pdb_dir, profile_dir)
    return rule

