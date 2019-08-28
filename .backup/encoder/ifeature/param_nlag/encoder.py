from typing import List, Callable, Optional, Dict, Any, Tuple, Iterable, Union

import numpy as np
import pandas as pd

from iFeature.codes.Moran import Moran
from iFeature.codes.Geary import Geary
from iFeature.codes.NMBroto import NMBroto
from iFeature.codes.SOCNumber import SOCNumber
from iFeature.codes.QSOrder import QSOrder

from encoder.encoder import BaseEncoder


class _ParamNLagEncoder(BaseEncoder):

    def _encode(self) -> Union[pd.DataFrame, Iterable[pd.DataFrame]]:
        res = self.func(self.in_data[0], nlag=self.nlag, **dict({"order": None}))
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name("ifeature", type(self).__name__.lower(), **{"nlag": self.nlag})
        )

    def encode(self) -> Iterable[pd.DataFrame]:
        if type(self.nlag) == range:
            for g in self.nlag:
                self.nlag = g
                yield self._encode()
        else:
            yield self._encode()

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]], int, Optional[Dict[str, Any]]], Tuple[List[str], List[float]]],
            nlag: int):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        if not self._sequences_sufficent_length(self.sequences, nlag, min_len=lambda _nlag: _nlag+1):
            raise ValueError(func.__name__ + ": all sequences should be greater than (nlag+1).")
        self.func = func
        self.nlag = nlag # default nlag: inspect.getcallargs(Moran, fastas=None)["nlag"]


class MoranEncoder(_ParamNLagEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, nlag: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, Moran, nlag)


class GearyEncoder(_ParamNLagEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, nlag: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, Geary, nlag)


class NMBrotoEncoder(_ParamNLagEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, nlag: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, NMBroto, nlag)


class SOCNumberEncoder(_ParamNLagEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, nlag: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, SOCNumber, nlag)


class QSOrderEncoder(_ParamNLagEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, nlag: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, QSOrder, nlag)
