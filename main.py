#!/usr/bin/env python

from snakemake import shell
from snakemake.io import expand

import nodes.utils as utils
import nodes.encodings as encodings

import sys

from peptidereactor.workflow_executer \
    import WorkflowExecuter, WorkflowSetter

CORES = 12
DATASETS = ["hiv_protease"]

with WorkflowSetter(cores=CORES) as w:

    w.add(utils.map_sequence_names.rule(
        fasta_in="data/{dataset}/seqs.fasta",
        fasta_out="data/{dataset}/seqs_mapped.fasta",
        maps_out="data/{dataset}/misc/mapped_sequence_names.yaml"))

    w.add(utils.tertiary_structure_search.rule(
        fasta_in="data/{dataset}/seqs_mapped.fasta", classes_in="data/{dataset}/classes.txt",
        fasta_sec_out="data/{dataset}/seqs_sec.fasta", classes_sec_out="data/{dataset}/classes_sec.txt",
        fasta_ter_out="data/{dataset}/seqs_ter.fasta", classes_ter_out="data/{dataset}/classes_ter.txt",
        pdb_dir="data/{dataset}/pdb/", profile_dir="data/{dataset}/profile/"))

    w.add(utils.multiple_sequence_alignment.rule(
        fastas_in=["data/{dataset}/seqs_mapped.fasta", "data/{dataset}/seqs_sec.fasta",
                   "data/{dataset}/seqs_ter.fasta"],
        fastas_out=["data/{dataset}/seqs_msa.fasta", "data/{dataset}/seqs_msa_sec.fasta",
                    "data/{dataset}/seqs_msa_ter.fasta"]))

    seqb = encodings.sequence_based.Rule()
    w.add(seqb.rule(
        fasta_in="data/{dataset}/seqs_mapped.fasta", fasta_msa_in="data/{dataset}/seqs_msa.fasta",
        classes_in="data/{dataset}/classes.txt", path_to_config="config.yaml",
        misc_dir="data/{dataset}/misc/", csv_dir="data/{dataset}/csv/"))

    strb = encodings.structure_based.Rule()
    w.add(strb.rule(
        include=["asa", "distance_distribution"],
        fasta_sec_in="data/{dataset}/seqs_sec.fasta", fasta_msa_sec_in="data/{dataset}/seqs_msa_sec.fasta",
        classes_sec_in="data/{dataset}/classes_sec.txt", fasta_ter_in="data/{dataset}/seqs_ter.fasta",
        classes_ter_in="data/{dataset}/classes_ter.txt", path_to_config="config.yaml", pdb_dir="data/{dataset}/pdb/",
        profile_dir="data/{dataset}/profile/", csv_dir="data/{dataset}/csv/"))

target = \
    expand(strb.target_csvs, dataset=DATASETS) + \
    expand(seqb.target_csvs, dataset=DATASETS)

target += expand(["data/{dataset}/seqs_msa.fasta", "data/{dataset}/seqs_msa_sec.fasta", "data/{dataset}/seqs_msa_ter.fasta"], dataset=DATASETS)
target += expand(["data/{dataset}/seqs_sec.fasta", "data/{dataset}/classes_sec.txt", "data/{dataset}/seqs_ter.fasta", "data/{dataset}/classes_ter.txt"], dataset=DATASETS)

with WorkflowExecuter(dict(), dict(out=target), "peptidereactor.yaml", cores=CORES) as e:
    shell(f"""./peptidereactor/run_pipeline -s peptidereactor.smk --configfile peptidereactor.yaml {" ".join(sys.argv[1:])}""")