from snakemake.io import expand

import yaml

import nodes.utils as utils
import nodes.encodings as encodings


class Rule:

    _MISC_ENCDOIGNS = \
        ["zscale", "dpc", "tpc", "gtpc", "gdpc", "gaac", "egaac", "dde",
         "ctdt", "ctdd", "ctdc", "blosum62", "binary", "aac", "blomap",
         "ctriad"]

    _PARAM_BASED_ENCODINGS = \
        ["ksctriad", "moran", "nmbroto", "geary", "qsorder", "socnumber",
         "eaac", "cksaagp", "cksaap", "apaac", "paac", "aaindex", "waac",
         "flgc", "fldpc", "fft", "cgr", "distance_frequency"]

    _NGRAM_BASED_ENCODINGS = \
        ["ngram_a2", "ngram_a3", "ngram_e2", "ngram_e3", "ngram_s2", "ngram_s3"]

    _PSEKRAAC_BASED_ENCODINGS = \
        [f"psekraac_type{t}" for t in
         [1, 2, 4, 5, ] + list(range(7, 17)) +
         ["3A", "3B", "6A", "6B", "6C"]]

    def _encodings_for_length_calculation(self, encodings=None):
        all = self._PARAM_BASED_ENCODINGS
        if encodings == None:
            return all
        else:
            return list(set(encodings).intersection(all))

    def _encodings_for_dim_calculation(self, encodings=None):
        all = self._NGRAM_BASED_ENCODINGS
        if encodings == None:
            return all
        else:
            return list(set(encodings).intersection(all))

    def _expand(self, src_dir, src, **wildcards):
        tmp = expand(src, **wildcards)
        return [f"{src_dir}{p}" for p in tmp]

    def rule(self, fasta_in, fasta_msa_in, classes_in, path_to_config, misc_dir, csv_dir, exclude=None, include=None):

        target_encodings = \
            self._MISC_ENCDOIGNS + \
            self._PARAM_BASED_ENCODINGS + \
            self._PSEKRAAC_BASED_ENCODINGS + \
            self._NGRAM_BASED_ENCODINGS

        if (exclude, type(include)) == (None, list):
            target_encodings = include
        elif (type(exclude), include) == (list, None):
            target_encodings = list(set(target_encodings).difference(exclude))
        elif (type(exclude), type(include)) == (list, list):
            print("Either param 'include' or 'exclude' must be greater zero. Ignoring the latter.")
            target_encodings = include

        with open(path_to_config) as f:
            config = yaml.safe_load(f)

        enc_len_calc, enc_dim_calc = \
            list(set(target_encodings).intersection(self._PARAM_BASED_ENCODINGS)), \
            list(set(target_encodings).intersection(self._NGRAM_BASED_ENCODINGS))

        rule = ""

        ### parameter computation

        if len(enc_len_calc) > 0:
            encodings_len_calc = self._encodings_for_length_calculation(enc_len_calc)
            rule += utils.window_length.rule(
                fasta_in, self._expand(misc_dir, "{encoding}.yaml", encoding=encodings_len_calc))

        if len(enc_dim_calc) > 0:
            encodings_dim_calc = self._encodings_for_dim_calculation(enc_dim_calc)
            rule += utils.dim_size.rule(
                fasta_in, self._expand(misc_dir, "{encoding}.yaml", encoding=encodings_dim_calc))

        ### misc encodings

        if "aac" in target_encodings \
                or "waac" in target_encodings \
                or "flgc" in target_encodings:
            aac_out = f"{csv_dir}aac.csv"
            rule += encodings.aac.rule(fasta_in, classes_in, aac_out)
            self.target_csvs += [aac_out]

        if "binary" in target_encodings:
            binary_out = f"{csv_dir}binary.csv"
            rule += encodings.binary.rule(fasta_in, classes_in, binary_out)
            self.target_csvs += [binary_out]

        if "blomap" in target_encodings:
            blomap_out = f"{csv_dir}blomap.csv"
            rule += encodings.blomap.rule(fasta_in, classes_in, blomap_out)
            self.target_csvs += [blomap_out]

        if "blosum62" in target_encodings:
            blosum62_out = f"{csv_dir}blosum62.csv"
            rule += encodings.blosum62.rule(fasta_in, classes_in, blosum62_out)
            self.target_csvs += [blosum62_out]

        if "ctdc" in target_encodings:
            ctdc_out = f"{csv_dir}ctdc.csv"
            rule += encodings.ctdc.rule(fasta_in, classes_in, ctdc_out)
            self.target_csvs += [ctdc_out]

        if "ctdd" in target_encodings:
            ctdd_out = f"{csv_dir}ctdd.csv"
            rule += encodings.ctdd.rule(fasta_in, classes_in, ctdd_out)
            self.target_csvs += [ctdd_out]

        if "ctdt" in target_encodings:
            ctdt_out = f"{csv_dir}ctdt.csv"
            rule += encodings.ctdt.rule(fasta_in, classes_in, ctdt_out)
            self.target_csvs += [ctdt_out]

        if "ctriad" in target_encodings:
            ctriad_out = f"{csv_dir}ctriad.csv"
            rule += encodings.ctriad.rule(fasta_in, classes_in, ctriad_out)
            self.target_csvs += [ctriad_out]

        if "dde" in target_encodings:
            dde_out = f"{csv_dir}dde.csv"
            rule += encodings.dde.rule(fasta_in, classes_in, dde_out)
            self.target_csvs += [dde_out]

        if "dpc" in target_encodings \
                or "fldpc" in target_encodings \
                or "ngram_a2" in target_encodings:
            dpc_out = f"{csv_dir}dpc.csv"
            rule += encodings.dpc.rule(fasta_in, classes_in, dpc_out)
            self.target_csvs += [dpc_out]

        if "egaac" in target_encodings:
            egaac_out = f"{csv_dir}egaac.csv"
            rule += encodings.egaac.rule(fasta_in, classes_in, egaac_out)
            self.target_csvs += [egaac_out]

        if "gaac" in target_encodings:
            gaac_out = f"{csv_dir}gaac.csv"
            rule += encodings.gaac.rule(fasta_in, classes_in, gaac_out)
            self.target_csvs += [gaac_out]

        if "gdpc" in target_encodings:
            gdpc_out = f"{csv_dir}gdpc.csv"
            rule += encodings.gdpc.rule(fasta_in, classes_in, gdpc_out)
            self.target_csvs += [gdpc_out]

        if "gtpc" in target_encodings:
            gtpc_out = f"{csv_dir}gtpc.csv"
            rule += encodings.gtpc.rule(fasta_in, classes_in, gtpc_out)
            self.target_csvs += [gtpc_out]

        if "tpc" in target_encodings or "ngram_a3" in target_encodings:
            tpc_out = f"{csv_dir}tpc.csv"
            rule = encodings.tpc.rule(fasta_in, classes_in, tpc_out)
            self.target_csvs += [tpc_out]

        if "zscale" in target_encodings:
            zscale_out = f"{csv_dir}zscale.csv"
            rule += encodings.zscale.rule(fasta_in, classes_in, zscale_out)
            self.target_csvs += [zscale_out]

        ### parameter-based encodings

        if "aaindex" in target_encodings or "fft" in target_encodings:
            aaindex_out = self._expand(csv_dir, "aaindex/aaindex_{aaindex}.csv", aaindex=config["aaindex"])
            rule += encodings.aaindex.rule(fasta_in, classes_in, aaindex_out)
            self.target_csvs += aaindex_out

        if "fft" in target_encodings:
            fft_out = self._expand(csv_dir, "fft/fft_{aaindex}.csv", aaindex=config["aaindex"])
            rule += encodings.fft.rule(aaindex_out, fft_out)
            self.target_csvs += fft_out

        if "waac" in target_encodings:
            waac_out = self._expand(csv_dir, "waac/waac_{aaindex}.csv", aaindex=config["aaindex"])
            rule += encodings.waac.rule(aac_out, waac_out)
            self.target_csvs += waac_out

        if "flgc" in target_encodings:
            flgc_out = self._expand(csv_dir, "flgc/flgc_{aaindex}.csv", aaindex=config["aaindex"])
            rule += encodings.flgc.rule(aac_out, flgc_out)
            self.target_csvs += flgc_out

        if "fldpc" in target_encodings:
            fldpc_out = self._expand(csv_dir, "fldpc/fldpc_{aaindex}.csv", aaindex=config["aaindex"])
            rule += encodings.fldpc.rule(dpc_out, fldpc_out)
            self.target_csvs += fldpc_out

        if "cgr" in target_encodings:
            cgr_out = self._expand(csv_dir, "cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
                                   resolution=config["cgr"]["resolution"], sfactor=config["cgr"]["sfactor"])
            rule += encodings.cgr.rule(fasta_in, classes_in, cgr_out)
            self.target_csvs += cgr_out

        if "distance_frequency" in target_encodings:
            distance_frequency_out = self._expand(csv_dir,
                                                  "distance_frequency/dist_freq_dn_{nterminal}_dc_{cterminal}.csv",
                                                  nterminal=config["distance_frequency"]["nterminal"],
                                                  cterminal=config["distance_frequency"]["cterminal"])
            rule += encodings.distance_frequency.rule(fasta_in, classes_in, distance_frequency_out)
            self.target_csvs += distance_frequency_out

        if "cksaagp" in target_encodings:
            cksaagp_out = self._expand(csv_dir, "cksaagp/cksaagp_gap_{gap_val}.csv", gap_val=range(*config["cksaagp"]))
            rule += encodings.cksaagp.rule(fasta_in, classes_in, f"{misc_dir}cksaagp.yaml", cksaagp_out)
            self.target_csvs += cksaagp_out

        ### psekraac-based encodings

        for t in [i.split("_")[1] for i in target_encodings if "psekraac" in i]:
            psekraac_out = \
                self._expand(csv_dir,
                             f"psekraac_{t}/{t.replace('ype', '')}_st-{{sub_val}}_rt-{{raac_val}}_ktu-{{ktuple_val}}_la-{{lambda_val}}.csv",
                             sub_val=config[f"psekraac_{t}"]["subtypes"],
                             raac_val=config[f"psekraac_{t}"]["raactypes"],
                             ktuple_val=config[f"psekraac_{t}"]["ktuples"],
                             lambda_val=config[f"psekraac_{t}"]["glambdas"])
            rule += eval(f"encodings.psekraac_{t}.rule(fasta_in, classes_in, psekraac_out)")
            self.target_csvs += psekraac_out

        ### ngram-based encodings

        for t in [i.split("_")[1] for i in target_encodings if "ngram" in i]:
            if t == "a2":
                ngram_a2_out = self._expand(csv_dir, f"ngram_{t}/ngram_{t}_{{dim}}.csv", dim=config[f"ngram_{t}"])
                ngram_a2_lsv_out = self._expand(misc_dir, f"ngram_{t}/ngram_{t}_lsv_{{dim}}.csv",
                                                dim=config[f"ngram_{t}"])
                ngram_a2_sv_out = self._expand(misc_dir, f"ngram_{t}/ngram_{t}_sv_{{dim}}.csv",
                                               dim=config[f"ngram_{t}"])
                rule += \
                    encodings.ngram.rule(t, csv_in=dpc_out, length_in=f"{misc_dir}ngram_{t}.yaml",
                                         ngram_out=ngram_a2_out, ngram_lsv_out=ngram_a2_lsv_out,
                                         ngram_sv_out=ngram_a2_sv_out)
                self.target_csvs += ngram_a2_out

            elif t == "a3":
                ngram_a3_out = self._expand(csv_dir, f"ngram_{t}/ngram_{t}_{{dim}}.csv", dim=config[f"ngram_{t}"])
                ngram_a3_lsv_out = self._expand(misc_dir, f"ngram_{t}/ngram_{t}_lsv_{{dim}}.csv",
                                                dim=config[f"ngram_{t}"])
                ngram_a3_sv_out = self._expand(misc_dir, f"ngram_{t}/ngram_{t}_sv_{{dim}}.csv",
                                               dim=config[f"ngram_{t}"])
                rule += \
                    encodings.ngram.rule(t, csv_in=tpc_out, length_in=f"{misc_dir}ngram_{t}.yaml",
                                         ngram_out=ngram_a3_out, ngram_lsv_out=ngram_a3_lsv_out,
                                         ngram_sv_out=ngram_a3_sv_out)
                self.target_csvs += ngram_a3_out

            else:
                ngram_out = self._expand(csv_dir, f"ngram_{t}/ngram_{t}_{{dim}}.csv", dim=config[f"ngram_{t}"])
                ngram_lsv_out = self._expand(misc_dir, f"ngram_{t}/ngram_{t}_lsv_{{dim}}.csv",
                                             dim=config[f"ngram_{t}"])
                ngram_sv_out = self._expand(misc_dir, f"ngram_{t}/ngram_{t}_sv_{{dim}}.csv",
                                            dim=config[f"ngram_{t}"])
                rule += \
                    encodings.ngram.rule(t, fasta_in=fasta_in, classes_in=classes_in,
                                         length_in=f"{misc_dir}ngram_{t}.yaml", ngram_out=ngram_out,
                                         ngram_lsv_out=ngram_lsv_out, ngram_sv_out=ngram_sv_out)

        return rule

    def __init__(self):
        self.target_csvs = []
