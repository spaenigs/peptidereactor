from snakemake import shell
from snakemake.utils import Namedlist
import yaml
import secrets
import os


class SnakemakeConfig:

    @staticmethod
    def get_cores(cores):
        return os.cpu_count() if cores == -1 else cores

    def _set_config(self):
        for obj in [self.input_files, self.output_files]:
            for k, v in obj.items():
                if type(v) == Namedlist:
                    obj[k] = list(v)
        return {**self.input_files, **self.output_files,
                **{"token": secrets.token_hex(4)}}

    def dump(self, path_to_configfile):
        with open(path_to_configfile, mode="w") as stream:
            yaml.safe_dump(self.config, stream)

    def __init__(self, input_files, output_files):
        self.input_files = input_files
        self.output_files = output_files
        self.config = self._set_config()


class MetaWorkflow(list):

    def register(self, wf):
        self.append(wf)

    def dryrun(self):
        for wf in self:
            shell(wf + " -nr")

