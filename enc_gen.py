#!/usr/bin/env python

from snakemake import shell
from snakemake.io import expand

import nodes.utils as utils
import nodes.encodings as encodings
import nodes.filter as dataset_filter
import nodes.benchmark as benchmark

import sys
import secrets

from peptidereactor.workflow_executer \
    import WorkflowExecuter, WorkflowSetter

TOKEN = secrets.token_hex(6)

CORES = 32
DATASETS = ["hiv_protease"]

with WorkflowSetter(cores=CORES, benchmark_dir="data/{dataset}/misc/benchmark/") as w:

    w.add(utils.check_dataset.rule(
        fasta_in="data/{dataset}/seqs.fasta", classes_in="data/{dataset}/classes.txt",
        report_out="data/{dataset}/misc/seqs_report.txt", benchmark_dir=w.benchmark_dir))

    w.add(utils.map_sequence_names.rule(
        fasta_in="data/{dataset}/seqs.fasta", fasta_out="data/{dataset}/seqs_mapped.fasta",
        maps_out="data/{dataset}/misc/mapped_sequence_names.yaml", benchmark_dir=w.benchmark_dir))

    w.add(utils.tertiary_structure_search.rule(
        fasta_in="data/{dataset}/seqs_mapped.fasta", classes_in="data/{dataset}/classes.txt",
        fasta_sec_out="data/{dataset}/seqs_sec.fasta", classes_sec_out="data/{dataset}/classes_sec.txt",
        fasta_ter_out="data/{dataset}/seqs_ter.fasta", classes_ter_out="data/{dataset}/classes_ter.txt",
        pdb_dir="data/{dataset}/pdb/", profile_dir="data/{dataset}/profile/", benchmark_dir=w.benchmark_dir))

    w.add(utils.multiple_sequence_alignment.rule(
        fastas_in=["data/{dataset}/seqs_mapped.fasta", "data/{dataset}/seqs_sec.fasta",
                   "data/{dataset}/seqs_ter.fasta"],
        fastas_out=["data/{dataset}/seqs_msa.fasta", "data/{dataset}/seqs_msa_sec.fasta",
                    "data/{dataset}/seqs_msa_ter.fasta"],
        benchmark_dir=w.benchmark_dir))

    seqb = encodings.sequence_based.Rule()
    w.add(seqb.rule(
        fasta_in="data/{dataset}/seqs_mapped.fasta", fasta_msa_in="data/{dataset}/seqs_msa.fasta",
        classes_in="data/{dataset}/classes.txt", path_to_config="config.yaml",
        misc_dir="data/{dataset}/misc/", csv_dir="data/{dataset}/csv/", benchmark_dir=w.benchmark_dir))

    strb = encodings.structure_based.Rule()
    w.add(strb.rule(
        fasta_sec_in="data/{dataset}/seqs_sec.fasta", fasta_msa_sec_in="data/{dataset}/seqs_msa_sec.fasta",
        classes_sec_in="data/{dataset}/classes_sec.txt", fasta_ter_in="data/{dataset}/seqs_ter.fasta",
        classes_ter_in="data/{dataset}/classes_ter.txt", path_to_config="config.yaml", pdb_dir="data/{dataset}/pdb/",
        profile_dir="data/{dataset}/profile/", csv_dir="data/{dataset}/csv/", benchmark_dir=w.benchmark_dir))

    w.add(utils.collect_encodings.rule(
        csv_seq_in=seqb.target_csvs, csv_str_in=strb.target_csvs,
        csv_seq_out=f"data/temp/{TOKEN}/{{dataset}}/csv/original/sequence_based/",
        csv_str_out=f"data/temp/{TOKEN}/{{dataset}}/csv/original/structure_based/"))

    sequence_based_encodings_dir, structure_based_encodings_dir = \
        "data/{dataset}/csv/sequence_based/", "data/{dataset}/csv/structure_based/"

    w.add(dataset_filter.non_empty.rule(
        csv_in=f"data/temp/{TOKEN}/{{dataset}}/csv/original/sequence_based/",
        csv_out=f"data/temp/{TOKEN}/{{dataset}}/csv/sequence_based/non_empty/", benchmark_dir=w.benchmark_dir))

    w.add(dataset_filter.curse_of_dim.rule(
        csv_in=f"data/temp/{TOKEN}/{{dataset}}/csv/sequence_based/non_empty/",
        csv_out=f"data/temp/{TOKEN}/{{dataset}}/csv/sequence_based/curse_of_dim/", benchmark_dir=w.benchmark_dir))

    w.add(dataset_filter.aaindex.rule(
        csv_in=f"data/temp/{TOKEN}/{{dataset}}/csv/sequence_based/curse_of_dim/",
        csv_out=f"data/temp/{TOKEN}/{{dataset}}/csv/sequence_based/aaindex/", benchmark_dir=w.benchmark_dir))

    w.add(dataset_filter.psekraac.rule(
        csv_in=f"data/temp/{TOKEN}/{{dataset}}/csv/sequence_based/aaindex/",
        csv_out=sequence_based_encodings_dir, benchmark_dir=w.benchmark_dir))

    w.add(dataset_filter.non_empty.rule(
        csv_in=f"data/temp/{TOKEN}/{{dataset}}/csv/original/structure_based/",
        csv_out=f"data/temp/{TOKEN}/{{dataset}}/csv/structure_based/non_empty/", benchmark_dir=w.benchmark_dir))

    w.add(dataset_filter.curse_of_dim.rule(
        csv_in=f"data/temp/{TOKEN}/{{dataset}}/csv/structure_based/non_empty/",
        csv_out=structure_based_encodings_dir, benchmark_dir=w.benchmark_dir))

    w.benchmark_target = [sequence_based_encodings_dir, structure_based_encodings_dir]
    target = \
        expand("data/{dataset}/misc/seqs_report.txt", dataset=DATASETS) + \
        expand(w.benchmark_dir + "benchmark.csv", dataset=DATASETS)

with WorkflowExecuter(dict(), dict(out=target), "peptidereactor.yaml", cores=CORES) as e:
    main_cmd = "./peptidereactor/run_pipeline -s peptidereactor.smk --configfile peptidereactor.yaml"
    options = f""" --cores {CORES} --keep-going {" ".join(sys.argv[1:])}"""
    shell(main_cmd + options)
