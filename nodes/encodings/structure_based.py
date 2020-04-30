from snakemake.io import expand

import yaml

import nodes.encodings as encodings


class Rule:

    # exclude for now disorder and pssm (see #15 and #18)
    _ENCODINGS = ["delaunay", "asa", "ta", "ssec", "sseb", "disorderb", "disorderc",
                  "qsar", "electrostatic_hull", "distance_distribution"]


    def _expand(self, src_dir, src, **wildcards):
        tmp = expand(src, **wildcards)
        return [f"{src_dir}{p}" for p in tmp]

    def rule(self,
             fasta_sec_in=None, fasta_msa_sec_in=None, classes_sec_in=None,
             fasta_ter_in=None, fasta_msa_ter_in=None, classes_ter_in=None,
             path_to_config=None, pdb_dir=None, profile_dir=None, csv_dir=None,
             exclude=None, include=None
             ):

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

        ### secondary-structure-based encodings

        if "asa" in target_encodings:
            asa_out = f"{csv_dir}asa.csv"
            rule += encodings.asa.rule(fasta_sec_in, classes_sec_in, profile_dir, asa_out)
            self.target_csvs += [asa_out]

        if "disorder" in target_encodings:
            disorder_out = f"{csv_dir}disorder.csv"
            rule += encodings.disorder.rule(fasta_sec_in, classes_sec_in, profile_dir, disorder_out)
            self.target_csvs += [disorder_out]

        if "disorderb" in target_encodings:
            disorderb_out = f"{csv_dir}disorderb.csv"
            rule += encodings.disorderb.rule(fasta_msa_sec_in, classes_sec_in, profile_dir, disorderb_out)
            self.target_csvs += [disorderb_out]

        if "disorderc" in target_encodings:
            disorderc_out = f"{csv_dir}disorderc.csv"
            rule += encodings.disorderc.rule(fasta_sec_in, classes_sec_in, profile_dir, disorderc_out)
            self.target_csvs += [disorderc_out]

        if "pssm" in target_encodings:
            pssm_out = f"{csv_dir}pssm.csv"
            rule += encodings.pssm.rule(fasta_sec_in, classes_sec_in, profile_dir, pssm_out)
            self.target_csvs += [pssm_out]

        if "sseb" in target_encodings:
            sseb_out = f"{csv_dir}sseb.csv"
            rule += encodings.sseb.rule(fasta_msa_sec_in, classes_sec_in, profile_dir, sseb_out)
            self.target_csvs += [sseb_out]

        if "ssec" in target_encodings:
            ssec_out = f"{csv_dir}ssec.csv"
            rule += encodings.ssec.rule(fasta_sec_in, classes_sec_in, profile_dir, ssec_out)
            self.target_csvs += [ssec_out]

        if "ta" in target_encodings:
            ta_out = f"{csv_dir}ta.csv"
            rule += encodings.ta.rule(fasta_sec_in, classes_sec_in, profile_dir, ta_out)
            self.target_csvs += [ta_out]

        ### tertiary-structure-based encodings

        if "delaunay" in target_encodings:
            delaunay_out = self._expand(csv_dir, "delaunay/delaunay_{algorithm}.csv", algorithm=config["delaunay"])
            rule += encodings.delaunay.rule(fasta_ter_in, classes_ter_in, pdb_dir, delaunay_out)
            self.target_csvs += delaunay_out

        if "distance_distribution" in target_encodings:
            distance_distribution_out = f"{csv_dir}distance_distribution.csv"
            rule += encodings.distance_distribution.rule(fasta_ter_in, classes_ter_in, pdb_dir, distance_distribution_out)
            self.target_csvs += [distance_distribution_out]

        if "electrostatic_hull" in target_encodings:
            electrostatic_hull_out = \
                self._expand(csv_dir, "electrostatic_hull/electrostatic_hull_{distance}.csv",
                             distance=config["electrostatic_hull"])
            rule += encodings.electrostatic_hull.rule(fasta_ter_in, classes_ter_in, pdb_dir, electrostatic_hull_out)
            self.target_csvs += electrostatic_hull_out

        if "qsar" in target_encodings:
            qsar_out = f"{csv_dir}qsar.csv"
            rule += encodings.qsar.rule(fasta_ter_in, classes_ter_in, pdb_dir, qsar_out)
            self.target_csvs += [qsar_out]

        return rule

    def __init__(self):
        self.target_csvs = []