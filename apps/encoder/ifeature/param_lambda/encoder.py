from typing import List, Callable, Optional, Dict, Any, Tuple, Iterable, Union

import numpy as np
import pandas as pd

from iFeature.codes.PAAC import PAAC
from iFeature.codes.APAAC import APAAC

from encoder.encoder import BaseEncoder


class _ParamLambdaEncoder(BaseEncoder):

    def _encode(self) -> Union[pd.DataFrame, Iterable[pd.DataFrame]]:
        res = self.func(self.in_data[0], lambdaValue=self.lambdaValue, **dict({"order": None}))
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name("ifeature", type(self).__name__.lower(), **{"lambda": self.lambdaValue})
        )

    def encode(self) -> Iterable[pd.DataFrame]:
        if type(self.lambdaValue) == range:
            for g in self.lambdaValue:
                self.lambdaValue = g
                yield self._encode()
        else:
            yield self._encode()

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]], int, Optional[Dict[str, Any]]], Tuple[List[str], List[float]]],
            lambdaValue: int):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        if not self._sequences_sufficent_length(self.sequences, lambdaValue, min_len=lambda _lV: _lV + 1):
            raise ValueError(func.__name__ + ": all sequences should be greater than (lambdaValue+1).")
        self.func = func
        self.lambdaValue = lambdaValue


class PAACEncoder(_ParamLambdaEncoder):

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            lambdaValue: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, PAAC, lambdaValue)


class APAACEncoder(_ParamLambdaEncoder):

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            lambdaValue: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, APAAC, lambdaValue)
