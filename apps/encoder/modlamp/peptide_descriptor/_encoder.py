from typing import Iterable, Union, List, Any, Type, Callable

import pandas as pd
import functools as ft
import inspect

from modlamp.descriptors import PeptideDescriptor

from encoder.encoder import BaseEncoder


def get_scales() -> List[str]:
    return [
        'aasi', 'abhprk', 'argos', 'bulkiness', 'charge_phys', 'charge_acid', 'cougar', 'eisenberg', 'ez',
        'flexibility', 'grantham', 'gravy', 'hopp-woods', 'isaeci', 'janin', 'kytedoolittle', 'levitt_alpha',
        'mss', 'msw', 'pepcats', 'peparc', 'polarity', 'ppcali', 'refractivity', 't_scale', 'tm_tend', 'z3', 'z5'
    ]


def get_binary_scales() -> List[str]:
    return ['abhprk', 'pepcats', 'peparc']


def get_profile_scales() -> List[str]:
    return [
        'aasi', 'argos', 'bulkiness', 'charge_phys', 'charge_acid', 'eisenberg',
        'flexibility', 'gravy', 'hopp-woods', 'janin', 'kytedoolittle', 'levitt_alpha',
        'mss', 'polarity', 'refractivity', 'tm_tend'
    ]


class PeptideDescriptorEncoder(BaseEncoder):

    def _encode(
            self,
            peptide_desc_obj: ft.partial,  # [PeptideDescriptor]
            func: Callable[..., int],
            windows: List[int],
            scales: List[str],
            modalities: List[str],
            prof_types: List[str]) -> Iterable[pd.DataFrame]:
        for w in windows:
            for s in scales:
                for m in modalities:
                    for pt in prof_types:
                        desc: PeptideDescriptor = peptide_desc_obj(scalename=s)
                        kwargs = {}
                        for kw in inspect.getfullargspec(func).args[1:-1]:
                            if kw == "window":
                                kwargs["window"] = w
                            elif kw == "modality":
                                kwargs["modality"] = m
                            elif kw == "prof_type":
                                kwargs["prof_type"] = pt
                            else:
                                raise ValueError("Unknown keyword: " + kw)
                        eval("desc." + func.__name__ + "(**{})".format(kwargs))
                        cn = map(
                            lambda x: desc.scalename.upper() + ".WINDOW" + str(w) + "." + str(x + 1),
                            range(desc.descriptor.shape[1])
                        )
                        yield self.to_df(
                            x=desc.descriptor,
                            y=self.in_data[1],
                            col_names=cn,
                            row_names=self.sequences,
                            meta_name=self.get_meta_name(
                                package="modlamp",
                                encoder=type(self).__name__.lower(),
                                **{"window": w, "scale": s, "modality": m, "proftype": pt}
                            )
                        )

    @staticmethod
    def _make_iterables(value: Union[Any, List[Any]], _type: Union[Type[int], Type[str]]):
        return [value] if type(value) == _type else value

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            window: Union[int, Iterable[int]],
            scale: Union[str, Iterable[str]]):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        self.window = window
        self.scale = scale


class AutocorrEncoder(PeptideDescriptorEncoder):

    def encode(self) -> Iterable[pd.DataFrame]:
        tmp_window = self._make_iterables(self.window, int)
        tmp_scale = self._make_iterables(self.scale, str)
        desc_partial = ft.partial(PeptideDescriptor, seqs=self.sequences)
        return super()._encode(
            peptide_desc_obj=desc_partial,
            func=PeptideDescriptor.calculate_autocorr,
            windows=tmp_window,
            scales=tmp_scale,
            modalities=[""],
            prof_types=[""]
        )


class CrosscorrEncoder(PeptideDescriptorEncoder):

    def encode(self) -> Iterable[pd.DataFrame]:
        tmp_window = self._make_iterables(self.window, int)
        tmp_scale = self._make_iterables(self.scale, str)
        desc_partial = ft.partial(PeptideDescriptor, seqs=self.sequences)
        return super()._encode(
            peptide_desc_obj=desc_partial,
            func=PeptideDescriptor.calculate_crosscorr,
            windows=tmp_window,
            scales=tmp_scale,
            modalities=[""],
            prof_types=[""]
        )


class MomentEncoder(PeptideDescriptorEncoder):
    # TODO angle in which to calculate the moment. 100 for alpha helices, 180 for beta sheets.
    # connect with ss prediction from iFeature?
    pass


class GlobalEncoder(PeptideDescriptorEncoder):

    def encode(self) -> Iterable[pd.DataFrame]:
        tmp_window = self._make_iterables(self.window, int)
        tmp_scale = self._make_iterables(self.scale, str)
        tmp_modality = self._make_iterables(self.modality, str)
        desc_partial = ft.partial(PeptideDescriptor, seqs=self.sequences)
        return super()._encode(
            peptide_desc_obj=desc_partial,
            func=PeptideDescriptor.calculate_global,
            windows=tmp_window,
            scales=tmp_scale,
            modalities=tmp_modality,
            prof_types=[""]
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            window: Union[int, Iterable[int]],
            scale: Union[str, Iterable[str]],
            modality: Union[str, Iterable[str]]):
        super().__init__(in_data, cores, window, scale)
        self.modality = modality


class ProfileEncoder(PeptideDescriptorEncoder):

    def encode(self) -> Iterable[pd.DataFrame]:
        tmp_window = self._make_iterables(self.window, int)
        tmp_scale = self._make_iterables(self.scale, str)
        tmp_prof_type = self._make_iterables(self.prof_type, str)
        desc_partial = ft.partial(PeptideDescriptor, seqs=self.sequences)
        return super()._encode(
            peptide_desc_obj=desc_partial,
            func=PeptideDescriptor.calculate_profile,
            windows=tmp_window,
            scales=tmp_scale,
            modalities=[""],
            prof_types=tmp_prof_type
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            window: Union[int, Iterable[int]],
            scale: Union[str, Iterable[str]],
            prof_type: Union[str, Iterable[str]]):
        super().__init__(in_data, cores, window, scale)
        self.prof_type = prof_type


class ArcEncoder(PeptideDescriptorEncoder):

    def encode(self) -> Iterable[pd.DataFrame]:
        tmp_scale = self._make_iterables(self.scale, str)
        tmp_modality = self._make_iterables(self.modality, str)
        desc_partial = ft.partial(PeptideDescriptor, seqs=self.sequences)
        return super()._encode(
            peptide_desc_obj=desc_partial,
            func=PeptideDescriptor.calculate_global,
            windows=[],
            scales=tmp_scale,
            modalities=tmp_modality,
            prof_types=[""]
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            scale: Union[str, Iterable[str]],
            modality: Union[str, Iterable[str]]):
        if not all([bs in get_binary_scales() for bs in self._make_iterables(scale, str)]):
            raise ValueError("Error: scales for ArcEncoder must be binary only.")
        super().__init__(in_data, cores, 1, scale)
        self.modality = modality




