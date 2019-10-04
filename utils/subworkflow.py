from functools import partial

import yaml
from snakemake import Workflow
from snakemake.workflow import Subworkflow

from snakemake import snakemake as smk_func


def set_yaml(_configfile, _input, _output):
    with open(_configfile, mode="w") as stream:
        yaml.safe_dump({**_input, **_output}, stream)


GLOBAL_WORKDIR = snakemake.config["global_workdir"]
NAME = snakemake.params.subworkflow
SNAKEFILE = snakemake.params.snakefile
CONFIGFILE = snakemake.params.configfile
INPUT, OUTPUT = dict(snakemake.input), dict(snakemake.output)

set_yaml(CONFIGFILE, INPUT, OUTPUT)

workflow = Workflow(GLOBAL_WORKDIR)
sw = Subworkflow(workflow, name=NAME, snakefile=SNAKEFILE, workdir=GLOBAL_WORKDIR, configfile=CONFIGFILE)
workflow._subworkflows[sw.name] = sw
workflow.globals[sw.name] = sw.target


@workflow.rule(name='all')
@workflow.input(sw.target(OUTPUT.values()))
@workflow.norun()
@workflow.run
def __rule_all(input, output, params, wildcards, threads, resources, log, version, rule, conda_env, singularity_img,
               singularity_args, use_singularity, bench_record, jobid, is_shell, bench_iteration, shadow_dir):
    pass

# os.chdir(snakemake.config["global_workdir"])
# workflow.w
workflow.check()
workflow.execute(dryrun=False, updated_files=[], quiet=True, resources=dict(),
                 subsnakemake=partial(smk_func))
