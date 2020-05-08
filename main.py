#!/usr/bin/env python

from snakemake import shell
from snakemake.io import expand

import nodes.utils as utils
import nodes.encodings as encodings
import nodes.filter as dataset_filter
import nodes.benchmark as benchmark

import sys

from peptidereactor.workflow_executer \
    import WorkflowExecuter, WorkflowSetter

CORES = 12
DATASETS = ["hiv_protease"]

with WorkflowSetter(cores=CORES) as w:

    w.add(utils.map_sequence_names.rule(
        fasta_in="data/{dataset}/seqs.fasta", fasta_out="data/{dataset}/seqs_mapped.fasta",
        maps_out="data/{dataset}/misc/mapped_sequence_names.yaml", benchmark_dir="data/{dataset}/misc/"))

    w.add(utils.tertiary_structure_search.rule(
        fasta_in="data/{dataset}/seqs_mapped.fasta", classes_in="data/{dataset}/classes.txt",
        fasta_sec_out="data/{dataset}/seqs_sec.fasta", classes_sec_out="data/{dataset}/classes_sec.txt",
        fasta_ter_out="data/{dataset}/seqs_ter.fasta", classes_ter_out="data/{dataset}/classes_ter.txt",
        pdb_dir="data/{dataset}/pdb/", profile_dir="data/{dataset}/profile/", benchmark_dir="data/{dataset}/misc/"))

    w.add(utils.multiple_sequence_alignment.rule(
        fastas_in=["data/{dataset}/seqs_mapped.fasta", "data/{dataset}/seqs_sec.fasta",
                   "data/{dataset}/seqs_ter.fasta"],
        fastas_out=["data/{dataset}/seqs_msa.fasta", "data/{dataset}/seqs_msa_sec.fasta",
                    "data/{dataset}/seqs_msa_ter.fasta"],
        benchmark_dir="data/{dataset}/misc/"))

    seqb = encodings.sequence_based.Rule()
    w.add(seqb.rule(
        fasta_in="data/{dataset}/seqs_mapped.fasta", fasta_msa_in="data/{dataset}/seqs_msa.fasta",
        classes_in="data/{dataset}/classes.txt", path_to_config="config.yaml",
        misc_dir="data/{dataset}/misc/", csv_dir="data/{dataset}/csv/", benchmark_dir="data/{dataset}/misc/"))

    strb = encodings.structure_based.Rule()
    w.add(strb.rule(
        fasta_sec_in="data/{dataset}/seqs_sec.fasta", fasta_msa_sec_in="data/{dataset}/seqs_msa_sec.fasta",
        classes_sec_in="data/{dataset}/classes_sec.txt", fasta_ter_in="data/{dataset}/seqs_ter.fasta",
        classes_ter_in="data/{dataset}/classes_ter.txt", path_to_config="config.yaml", pdb_dir="data/{dataset}/pdb/",
        profile_dir="data/{dataset}/profile/", csv_dir="data/{dataset}/csv/", benchmark_dir="data/{dataset}/misc/"))

    sequence_based_encodings_dir, structure_based_encodings_dir = \
        "data/{dataset}/csv/sequence_based/", "data/{dataset}/csv/structure_based/"

    w.add(dataset_filter.non_empty.rule(
        csv_in=seqb.target_csvs,
        csv_out="data/temp/{dataset}/csv/sequence_based/non_empty/",
        benchmark_dir="data/{dataset}/misc/benchmark/"))

    w.add(dataset_filter.curse_of_dim.rule(
        csv_in="lambda wildcards: glob('data/temp/{dataset}/csv/sequence_based/non_empty/*.csv')",
        csv_out="data/temp/{dataset}/csv/sequence_based/curse_of_dim/",
        benchmark_dir="data/{dataset}/misc/benchmark/"))

    w.add(dataset_filter.aaindex.rule(
        csv_in="lambda wildcards: glob(f'data/temp/{wildcards.dataset}/csv/sequence_based/curse_of_dim/*.csv')",
        csv_out="data/temp/{dataset}/csv/sequence_based/aaindex/", benchmark_dir="data/{dataset}/misc/benchmark/"))

    w.add(dataset_filter.psekraac.rule(
        csv_in="lambda wildcards: glob(f'data/temp/{wildcards.dataset}/csv/sequence_based/aaindex/*.csv')",
        csv_out=sequence_based_encodings_dir, benchmark_dir="data/{dataset}/misc/benchmark/"))

    w.add(dataset_filter.non_empty.rule(
        csv_in=strb.target_csvs,
        csv_out=structure_based_encodings_dir,
        benchmark_dir="data/{dataset}/misc/"))

    w.add(benchmark.cross_validation.rule(
        csv_in=["data/hiv_protease/csv/zscale.csv"],
        csv_dir_out="data/{dataset}/benchmark/",
        benchmark_dir="data/{dataset}/misc/benchmark/"))

target = \
    expand(sequence_based_encodings_dir, dataset=DATASETS) + \
    expand(structure_based_encodings_dir, dataset=DATASETS)

target += expand("data/{dataset}/benchmark/", dataset=DATASETS)

# target = expand(["data/{dataset}/seqs_msa.fasta", "data/{dataset}/seqs_msa_sec.fasta", "data/{dataset}/seqs_msa_ter.fasta"], dataset=DATASETS)
# target += expand(["data/{dataset}/seqs_sec.fasta", "data/{dataset}/classes_sec.txt", "data/{dataset}/seqs_ter.fasta", "data/{dataset}/classes_ter.txt"], dataset=DATASETS)

with WorkflowExecuter(dict(), dict(out=target), "peptidereactor.yaml", cores=CORES) as e:
    shell(f"""./peptidereactor/run_pipeline -s peptidereactor.smk --configfile peptidereactor.yaml {" ".join(sys.argv[1:])}""")