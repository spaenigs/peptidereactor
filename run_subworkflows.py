# from: https://gist.github.com/marcelm/7e432275dfa5e762b4f8
from functools import partial

from snakemake.workflow import Workflow, Rules, Subworkflow
import snakemake.dag

from utils.subworkflow import SubworkflowEB

def register(sw):
    workflow._subworkflows[sw.name] = sw
    workflow.globals[sw.name] = sw.target
    return sw.target(sw.output.values())


workflow = Workflow(__file__)
snakemake.workflow.rules = Rules()
snakemake.workflow.config = dict()

subworkflows = [
    SubworkflowEB(workflow, name="disorder_swf", workdir=".",
                  snakefile="nodes/encodings/disorder/Snakefile", configfile="nodes/encodings/disorder/config.yaml",
                  fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
                  classes_in="data/neuropeptides_ds3/annotated_classes.txt",
                  csv_out="data/neuropeptides_ds3/csv/disorder.csv"),
    SubworkflowEB(workflow, name="aac_swf", workdir=".",
                  snakefile="nodes/encodings/aac/Snakefile", configfile="nodes/encodings/aac/config.yaml",
                  fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
                  classes_in="data/neuropeptides_ds3/annotated_classes.txt",
                  csv_out="data/neuropeptides_ds3/csv/aac.csv")]


@workflow.rule(name='all')
@workflow.input([register(subworkflow) for subworkflow in subworkflows])
@workflow.norun()
@workflow.run
def __rule_all(input, output, params, wildcards, threads, resources, log, version, rule, conda_env, singularity_img,
               singularity_args, use_singularity, bench_record, jobid, is_shell, bench_iteration, shadow_dir):
    pass


workflow.check()
workflow.execute(dryrun=False, updated_files=[], quiet=True, resources=dict(),
                 subsnakemake=partial(snakemake.snakemake))
