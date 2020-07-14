from snakemake.utils import Namedlist
from shutil import rmtree

import yaml
import secrets
import os


class WorkflowExecuter:

    @staticmethod
    def get_cores(cores):
        return os.cpu_count() if cores == -1 else cores

    def _set_config(self):
        for obj in [self.input_files, self.output_files]:
            for k, v in obj.items():
                if type(v) == Namedlist:
                    obj[k] = list(v)
        return {**self.input_files, **self.output_files,
                **{"token": self.token}}

    def write_configfile(self):
        with open(self.path_to_configfile, mode="w") as stream:
            yaml.safe_dump(self.config, stream)

    def remove_configfile(self):
        os.remove(self.path_to_configfile)

    def remove_temp_dir(self):
        path = f"data/temp/{self.token}/"
        if os.path.exists(path):
            rmtree(path)

    def snakemake(self):
        return

    def __enter__(self):
        self.write_configfile()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.remove_configfile()
        # self.remove_temp_dir()

    def __init__(self, input_files, output_files, path_to_configfile, cores=1, **kwargs):
        self.snakemake = f"snakemake --nolock --quiet -d $PWD --cores {cores}"
        self.token = secrets.token_hex(6)
        self.input_files = input_files
        self.output_files = output_files
        self.config = self._set_config()
        self.path_to_configfile = path_to_configfile
        if len(kwargs) > 0:
            for key, value in kwargs.items():
                self.config[key] = value


class WorkflowSetter:

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        scaffold = f"""\
from peptidereactor.workflow_executer import WorkflowExecuter

CORES = {self.cores}

rule all:
    input:
         config['{self.key_target_rule}']
"""

        for r in self.rule_definitions:
            scaffold += r

        with open(self.snakefile, mode="w") as f:
            f.write(scaffold)

    def add(self, rule_definition):
        self.rule_definitions += [rule_definition]

    def __init__(self, cores=1, key_target_rule="out", snakefile="peptidereactor.smk", benchmark_dir=None):
        self.rule_definitions = []
        self.cores = cores
        self.key_target_rule = key_target_rule
        self.snakefile = snakefile
        self.benchmark_dir = benchmark_dir
        self.benchmark_target = []
