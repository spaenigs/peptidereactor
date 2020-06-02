from snakemake.utils import Namedlist
from shutil import rmtree

import yaml
import secrets
import os
import textwrap


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


class MetaWorkflowExecuter(WorkflowExecuter):

    def __init__(self, input_files, output_files, path_to_configfile, cores=1, **kwargs):
        super().__init__(input_files, output_files, path_to_configfile, cores, **kwargs)
        self.snakemake = f"snakemake --nolock --quiet -d $PWD --config cores={cores} --latency-wait 60"


class WorkflowSetter:

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        scaffold = textwrap.dedent(
            f"""\
                from glob import glob
            
                from peptidereactor.workflow_executer import WorkflowExecuter
                
                CORES = {self.cores}
                
                rule all:
                    input:
                         config['{self.key_target_rule}']
            """)

        for r in self.rule_definitions:
            scaffold += r

        if self.benchmark_dir is not None:
            scaffold += textwrap.dedent(
                f"""\

                    rule collect_benchmark:
                        input:
                             {self.benchmark_target}
                        output:
                             "{self.benchmark_dir}benchmark.csv"
                        threads:
                             1000
                        run:
                             import re
                             import pandas as pd

                             df_res = pd.DataFrame()
                             for p in glob(f"data/{{wildcards.dataset}}/misc/benchmark/*/*.txt"):
                                 name = re.findall(".*/(.*)_\w+.txt", p)[0]
                                 df_tmp = pd.read_csv(p, sep="\t")
                                 df_tmp.index = [name]
                                 df_res = pd.concat([df_res, df_tmp])

                             df_res.to_csv(output[0])
                """)

        with open(self.snakefile, mode="w") as f:
            f.write(scaffold)

    def add(self, rule_definition):
        self.rule_definitions += [rule_definition]

    def __init__(self, cores=1, key_target_rule="out", snakefile="peptidereactor.smk", benchmark_dir=None):
        self.rule_definitions = []
        self.cores = cores
        self.key_target_rule = key_target_rule
        self.snakefile = snakefile
        self.benchmark_dir = \
            None if benchmark_dir is None else f"{benchmark_dir}{secrets.token_hex(3)}/"
        self.benchmark_target = []
