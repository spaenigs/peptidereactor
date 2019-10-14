from functools import partial

import yaml
import os
import secrets
from snakemake import Workflow
from snakemake.workflow import Subworkflow
from snakemake.io import Namedlist

from snakemake import snakemake as smk_func


def set_input_output(_configfile, _input, _output, token):
    for k, v in _output.items():
        if type(v) == Namedlist:
            _output[k] = list(v)
    with open(_configfile, mode="a") as stream:
        yaml.safe_dump({**_input, **_output, **{"token": token}}, stream)


GLOBAL_WORKDIR = snakemake.config["global_workdir"]
NAME = snakemake.params.subworkflow
SNAKEFILE = snakemake.params.snakefile
CONFIGFILE = snakemake.params.configfile
INPUT, OUTPUT = dict(snakemake.input), dict(snakemake.output)
RESOURCES = dict(snakemake.resources)
TOKEN = secrets.token_hex(4)

set_input_output(CONFIGFILE, INPUT, OUTPUT, TOKEN)

workflow = Workflow(GLOBAL_WORKDIR, default_resources=None)
workflow.included_stack.append("")

sw = Subworkflow(workflow,
                 name=NAME,
                 snakefile=SNAKEFILE,
                 workdir=GLOBAL_WORKDIR,
                 configfile=CONFIGFILE)

workflow._subworkflows[sw.name] = sw
workflow.globals[sw.name] = sw.target


@workflow.rule(name='all')
@workflow.input(sw.target(OUTPUT.values()))
@workflow.norun()
@workflow.run
def __rule_all(input, output, params, wildcards, threads, resources, log, version, rule, conda_env, singularity_img,
               singularity_args, use_singularity, bench_record, jobid, is_shell, bench_iteration, shadow_dir):
    pass


use_cores = os.cpu_count() if RESOURCES.get("cores", 1) == -1 else RESOURCES.get("cores", 1)

workflow.check()
success = workflow.execute(dryrun=False, updated_files=[], quiet=True, resources=dict(),
                           subsnakemake=partial(smk_func, cores=use_cores, quiet=False, printreason=True))

os.remove(CONFIGFILE)
os.removedirs(f"data/temp/{TOKEN}")

