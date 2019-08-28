from typing import List, Callable, Optional, Dict, Any, Tuple, Iterable, Union

import numpy as np
import pandas as pd

from iFeature.codes.CKSAAGP import CKSAAGP
from iFeature.codes.CKSAAP import CKSAAP
from iFeature.codes.KSCTriad import KSCTriad
from iFeature.codes.CTriad import CTriad

from encoder.encoder import BaseEncoder


class _ParamGapEncoder(BaseEncoder):

    def _encode(self) -> Union[pd.DataFrame, Iterable[pd.DataFrame]]:
        res = self.func(self.in_data[0], self.gap, **dict({"order": None}))
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name("ifeature", type(self).__name__.lower(), **{"gap": self.gap})
        )

    def encode(self) -> Iterable[pd.DataFrame]:
        if type(self.gap) == range:
            for g in self.gap:
                self.gap = g
                yield self._encode()
        else:
            yield self._encode()

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]], int, Optional[Dict[str, Any]]], Tuple[List[str], List[float]]],
            gap: int):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        self.func = func
        self.gap = gap


class CKSAAPEncoder(_ParamGapEncoder):

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            gap: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, CKSAAP, gap)
        if not self._sequences_sufficent_length(self.sequences, gap, min_len=lambda _gap: _gap+2):
            raise ValueError("CKSAAP: all sequences should be greater than (gap+2).")


class CKSAAGPEncoder(_ParamGapEncoder):

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            gap: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, CKSAAGP, gap)
        if not self._sequences_sufficent_length(self.sequences, gap, min_len=lambda _gap: _gap+2):
            raise ValueError("CKSAAGP: all sequences should be greater than (gap+2).")


class KSCTriadEncoder(_ParamGapEncoder):

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            gap: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, KSCTriad, gap)
        if not self._sequences_sufficent_length(self.sequences, gap, min_len=lambda _gap: 2*_gap+3):
            raise ValueError("All sequences should be greater than (2*gap+3).")
        self.gap = gap


class CTriadEncoder(_ParamGapEncoder):

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            gap: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, CTriad, gap)
        if not self._sequences_sufficent_length(self.sequences, gap, min_len=lambda _gap: _gap+3):
            raise ValueError("All sequences should be greater than (gap+3).")
