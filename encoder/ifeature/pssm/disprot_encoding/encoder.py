from typing import List

import subprocess as sp
import pandas as pd
import numpy as np
import functools as ft

from iFeature.codes.Disorder import Disorder
from iFeature.codes.DisorderB import DisorderB
from iFeature.codes.DisorderC import DisorderC

from encoder.encoder import BaseEncoder, BaseEncoderEqualLength
from encoder.ifeature.pssm.utils import PSSMUtils


class _DisorderEncoder(PSSMUtils):

    @staticmethod
    def generate_profile(in_data, profile_dir, cores, **kwargs):

        for record in in_data[0]:
            sp.run("{} {} {} {} {} {}".format(
                "/home/spaenigs/Apps/VSL2/runvsl2",
                record[1],
                profile_dir,
                profile_dir,
                record[0],
                "/home/spaenigs/Apps/VSL2"
            ), shell=True, check=True)


class DisorderEncoder(BaseEncoderEqualLength, _DisorderEncoder):

    def encode(self) -> pd.DataFrame:
        # self.generate_profile(
        #     in_data=self.in_data,
        #     profile_dir=self.profile_dir,
        #     cores=self.cores
        # )
        return super()._encode(ft.partial(Disorder, path=self.profile_dir))

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            profile_dir: str):
        super().__init__(in_data, cores, 1)
        PSSMUtils.check_pssm_profile_presence(profile_dir)
        self.profile_dir = profile_dir


class DisorderBEncoder(BaseEncoderEqualLength, _DisorderEncoder):

    def encode(self) -> pd.DataFrame:
        # self.generate_profile(
        #     in_data=self.in_data,
        #     profile_dir=self.profile_dir,
        #     cores=self.cores
        # )
        return super()._encode_binary_with_gaps(
            fastas_aligned=self.fastas_aligned,
            func=ft.partial(DisorderB, path=self.profile_dir),
            binary_len=2
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            profile_dir: str,
            run_msa: bool):
        super().__init__(in_data, cores, 1)
        PSSMUtils.check_pssm_profile_presence(profile_dir)
        if run_msa:
            mean_seq_len, self.fastas_aligned = self.run_muscle(in_data)
        else:
            self.fastas_aligned = self.in_data[0]
            in_data_tmp_seq = []
            for tup in self.in_data[0]:
                in_data_tmp_seq.append([tup[0], tup[1].replace("-", "")])
            self.in_data = (in_data_tmp_seq, self.in_data[1])
            mean_seq_len = len(self.fastas_aligned[0][1])
        self.interpolate_to = mean_seq_len * 2
        self.profile_dir = profile_dir


class DisorderCEncoder(_DisorderEncoder, BaseEncoder):

    def encode(self) -> pd.DataFrame:
        # self.generate_profile(
        #     in_data=self.in_data,
        #     profile_dir=self.profile_dir,
        #     cores=self.cores
        # )
        res = DisorderC(self.in_data[0], path=self.profile_dir)
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name(package="ifeature", encoder=type(self).__name__.lower())
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            profile_dir: str):
        super().__init__(in_data, cores)
        PSSMUtils.check_pssm_profile_presence(profile_dir)
        self.profile_dir = profile_dir
        self.sequences = [tup[1] for tup in self.in_data[0]]
