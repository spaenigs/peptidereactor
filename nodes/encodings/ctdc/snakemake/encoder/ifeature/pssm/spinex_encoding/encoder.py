from typing import List

import pandas as pd
import functools as ft
import subprocess as sp
from pathos.multiprocessing import ProcessingPool as Pool

from iFeature.codes.ASA import ASA
from iFeature.codes.TA import TA

from encoder.ifeature.pssm.utils import PSSMUtils
from encoder.encoder import BaseEncoderEqualLength


class _SpineXEncoder(PSSMUtils):

    @staticmethod
    def generate_profile(in_data, profile_dir, cores, **kwargs):

        for record in in_data[0]:
            sp.run("{} {} {} {}".format(
                "/home/spaenigs/Apps/spineXpublic/runspx",
                "/home/spaenigs/Apps/spineXpublic",
                profile_dir,
                record[0]
            ), shell=True, check=True)


class ASAEncoder(BaseEncoderEqualLength, _SpineXEncoder):

    def encode(self) -> pd.DataFrame:
        # self.generate_profile(
        #     in_data=self.in_data,
        #     profile_dir=self.profile_dir,
        #     cores=self.cores
        # )
        return super()._encode(ft.partial(ASA, path=self.profile_dir))

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            profile_dir: str):
        super().__init__(in_data, cores, 1)
        PSSMUtils.check_pssm_profile_presence(profile_dir)
        self.profile_dir = profile_dir


class TAEncoder(BaseEncoderEqualLength, _SpineXEncoder):

    def encode(self) -> pd.DataFrame:
        # self.generate_profile(
        #     in_data=self.in_data,
        #     profile_dir=self.profile_dir,
        #     cores=self.cores
        # )
        return super()._encode(ft.partial(TA, path=self.profile_dir))

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            profile_dir: str):
        super().__init__(in_data, cores, 1)
        PSSMUtils.check_pssm_profile_presence(profile_dir)
        self.profile_dir = profile_dir
