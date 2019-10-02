from snakemake.workflow import Subworkflow
import yaml


class SubworkflowEB(Subworkflow):

    def __set_yaml(self):
        with open(self.configfile, mode="w") as stream:
            yaml.safe_dump({**self.input, **self.output}, stream)

    def __init__(self, workflow, name, snakefile, workdir, configfile, **kwargs):
        super().__init__(workflow, name, snakefile, workdir, configfile)
        self.input, self.output = {}, {}
        for key, value in kwargs.items():
            if "in" in key:
                self.input[key] = value
            elif "out" in key:
                self.output[key] = value
            else:
                raise ValueError("Specify in- and output files with _in or _out, respectively.")
        self.__set_yaml()

