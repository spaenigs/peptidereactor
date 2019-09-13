import glob
import os
import re
from typing import List, Iterator

from Bio.Blast.Applications import NcbipsiblastCommandline

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import glob as gl
import numpy as np
from shutil import copyfile


class PSSMUtils:

    @staticmethod
    def save_fasta(in_data: (List[List[str]], List[int]), filename: str) -> None:

        raw_data = in_data[0]
        target = in_data[1]

        def _get_records(_raw_data: List[List[str]], _target: List[int]) -> Iterator[SeqRecord]:
            for tup, y, i in zip(_raw_data, _target, range(1, 100000000)):
                cl, seq = tup
                yield SeqRecord(
                    seq=Seq(seq),
                    id="Seq_" + str(i),
                    description=str(y),
                    name="Seq_" + str(i)
                )

        SeqIO.write(_get_records(raw_data, target), filename, "fasta")

    @staticmethod
    def read_fasta(filename: str) -> List[List[str]]:
        raw_data: List[List[str]] = []
        for record in SeqIO.parse(filename, "fasta"):
            raw_data.append([record.name, str(record.seq)])
        return raw_data

    @staticmethod
    def filter_non_hits(
            in_data: (List[List[str]], List[int]),
            profile_dir: str,
            file_pattern: str) -> (List[List[str]], List[int]):
        raw_data = in_data[0]
        target = in_data[1]
        # res = sorted(
        #     [tup for f in gl.glob(profile_dir + file_pattern) for tup in raw_data if tup[0] in f],
        #     key=lambda tup: int(tup[0].replace("Seq_", ""))
        # )
        # idx = np.unique(res, return_index=True, axis=0)[1]
        # return [res[i] for i in idx], [target[i] for i in idx]
        unique_names = set(
            [re.sub("\.\w{1,}", "", os.path.basename(os.path.splitext(p)[0])) for p in
             glob.glob(profile_dir + file_pattern)]
        )
        non_hits = set(map(lambda tup: tup[0], raw_data)).difference(unique_names)
        zipped = list(filter(lambda tup: tup[0][0] not in non_hits, zip(raw_data, target)))
        raw_data_filtered = list(map(lambda tup: tup[0], zipped))
        target_filtered = list(map(lambda tup: tup[1], zipped))
        return raw_data_filtered, target_filtered


    @staticmethod
    def generate_profile(
            in_data: (List[List[str]], List[int]),
            profile_dir: str,
            cores: int,
            **kwargs) -> None:

        for record in in_data[0]:
            query_string = '>' + record[0] + '\n' + str(record[1])
            matrix_file_name = "{}/{}.pssm".format(profile_dir, record[0])
            cline = NcbipsiblastCommandline(
                cmd="psiblast",
                db=kwargs["db"],
                inclusion_ethresh=0.001,
                num_iterations=3,
                num_threads=cores,
                out_pssm="{}/{}.asn.pssm".format(profile_dir, record[0]),
                out_ascii_pssm=matrix_file_name
            )
            o, e = cline(stdin=query_string)
            try:
                copyfile(
                    matrix_file_name,
                    matrix_file_name.replace(".pssm", ".mat")
                )
            except Exception as e:
                print("e" + str(e))
                print("Psiblast: no hits found. Copy result files not necessary.")

    @staticmethod
    def check_pssm_profile_presence(profile_dir: str):
        files = [f for f in gl.glob(profile_dir + "/*.pssm")] + \
            [f for f in gl.glob(profile_dir + "/*.asn.pssm")] + \
            [f for f in gl.glob(profile_dir + "/*.mat")]
        if len(files) == 0:
            raise FileNotFoundError("No PSSM profile found. Please run PSSMUtils.generate_profile first.")

