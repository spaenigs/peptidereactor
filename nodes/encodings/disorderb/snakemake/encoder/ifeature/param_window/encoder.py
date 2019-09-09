from typing import List, Callable, Optional, Dict, Any, Tuple, Iterable, Union

import functools as ft
import pandas as pd

from iFeature.codes.EAAC import EAAC
from iFeature.codes.EGAAC import EGAAC

from encoder.encoder import BaseEncoderEqualLength


class _ParamWindowEncoder(BaseEncoderEqualLength):

    def encode(self) -> Iterable[pd.DataFrame]:
        if type(self.window) == range:
            for g in self.window:
                self.window = g
                yield self._encode(ft.partial(self.func, window=self.window, **dict({"order": None})), **{"window": self.window})
        else:
            yield self._encode(ft.partial(self.func, window=self.window, **dict({"order": None})), **{"window": self.window})

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]], int, Optional[Dict[str, Any]]], Tuple[List[str], List[float]]],
            window: int):
        super().__init__(in_data, cores, 1)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        if not self._sequences_sufficent_length(self.sequences, window, min_len=lambda _window: _window):
            raise ValueError(func.__name__ + ": all sequences should be greater than sliding window.")
        self.func = func
        self.window = window


class EAACEncoder(_ParamWindowEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, window: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, EAAC, window)


class EGAACEncoder(_ParamWindowEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, window: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, EGAAC, window)
