from typing import List

import functools as ft
import pandas as pd

from iFeature.codes.PSSM import PSSM

from encoder.encoder import BaseEncoderEqualLength
from encoder.ifeature.pssm.utils import PSSMUtils


class PSSMEncoder(BaseEncoderEqualLength):

    def encode(self) -> pd.DataFrame:
        return super()._encode(ft.partial(PSSM, path=self.profile_dir))

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            profile_dir: str):
        # in_data_filtered = PSSMUtils.filter_non_hits(
        #     in_data=in_data,
        #     profile_dir=profile_dir,
        #     file_pattern="/*asn.pssm"
        # )
        super().__init__(in_data, cores, 20)
        PSSMUtils.check_pssm_profile_presence(profile_dir)
        self.profile_dir = profile_dir
