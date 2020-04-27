from snakemake.io import expand

import yaml

import nodes.encodings as encodings

class Rule:

    _ENCODINGS = ["delaunay"]


    def _expand(self, src_dir, src, **wildcards):
        tmp = expand(src, **wildcards)
        return [f"{src_dir}{p}" for p in tmp]

    def rule(self, fasta_in, fasta_msa_in, classes_in, path_to_config, pdb_dir, profile_dir, csv_dir, exclude=None, include=None):

        target_encodings = self._ENCODINGS

        if (exclude, type(include)) == (None, list):
            target_encodings = include
        elif (type(exclude), include) == (list, None):
            target_encodings = list(set(target_encodings).difference(exclude))
        elif (type(exclude), type(include)) == (list, list):
            print("Either param 'include' or 'exclude' must be greater zero. Ignoring the latter.")
            target_encodings = include

        with open(path_to_config) as f:
            config = yaml.safe_load(f)

        rule = ""

        if "delaunay" in target_encodings:
            delaunay_out = self._expand(csv_dir, "delaunay/delaunay_{algorithm}.csv", algorithm=config["delaunay"])
            rule += encodings.delaunay.rule(fasta_in, classes_in, pdb_dir, delaunay_out)
            self.target_csvs += delaunay_out

        return rule

    def __init__(self):
        self.target_csvs = []