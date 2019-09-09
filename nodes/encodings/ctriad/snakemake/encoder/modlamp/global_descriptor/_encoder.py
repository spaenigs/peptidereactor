from typing import List

from modlamp.descriptors import GlobalDescriptor

import pandas as pd
import numpy as np

from encoder.encoder import BaseEncoder


class GlobalDescriptorEncoder(BaseEncoder):

    def encode(self) -> pd.DataFrame:
        global_desc = GlobalDescriptor(self.sequences)
        global_desc.calculate_all()
        return self.to_df(
            x=np.array(global_desc.descriptor),
            y=self.in_data[1],
            col_names= global_desc.featurenames,
            row_names=self.sequences
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int):
        super().__init__(in_data, cores)
        self.sequences = list(map(lambda tup: tup[1], in_data[0]))


__all__ = ["GlobalDescriptorEncoder"]
