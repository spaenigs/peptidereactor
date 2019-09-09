from typing import List, Callable, Dict, Tuple, Union, Iterable

import pandas as pd
import numpy as np

import iFeature.PseKRAAC.type1 as t1
import iFeature.PseKRAAC.type2 as t2
import iFeature.PseKRAAC.type3A as t3A
import iFeature.PseKRAAC.type3B as t3B
import iFeature.PseKRAAC.type4 as t4
import iFeature.PseKRAAC.type5 as t5
import iFeature.PseKRAAC.type6A as t6A
import iFeature.PseKRAAC.type6B as t6B
import iFeature.PseKRAAC.type6C as t6C
import iFeature.PseKRAAC.type7 as t7
import iFeature.PseKRAAC.type8 as t8
import iFeature.PseKRAAC.type9 as t9
import iFeature.PseKRAAC.type10 as t10
import iFeature.PseKRAAC.type11 as t11
import iFeature.PseKRAAC.type12 as t12
import iFeature.PseKRAAC.type13 as t13
import iFeature.PseKRAAC.type14 as t14
import iFeature.PseKRAAC.type15 as t15
import iFeature.PseKRAAC.type16 as t16

from encoder.encoder import BaseEncoder


class _PseKRAACEncoder(BaseEncoder):

    def _encode(self, subtype: str, raactype: int, ktuple: int, glambda: int):
        res = self.func(
            self.in_data[0],
            subtype=subtype,
            raactype=raactype,
            ktuple=ktuple,
            glValue=glambda
        )
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name(
                package="ifeature",
                encoder=type(self).__name__.lower(),
                **{"subtype": subtype, "raactype": raactype, "ktuple": ktuple, "glValue": glambda}
            )
        )

    def _make_iterables(self):
        tmp_subtype = \
            [self.subtype] if type(self.subtype) == str else self.subtype
        tmp_raactype = \
            [self.raactype] if type(self.raactype) == int else self.raactype
        tmp_ktuple = \
            [self.ktuple] if type(self.ktuple) == int else self.ktuple
        tmp_glambda = \
            [self.glambda] if type(self.glambda) == int else self.glambda
        return (tmp_subtype, tmp_raactype, tmp_ktuple, tmp_glambda)

    def encode(self) -> Iterable[pd.DataFrame]:
        tmp_subtype, tmp_raactype, tmp_ktuple, tmp_glambda = self._make_iterables()
        for st in tmp_subtype:
            for rt in tmp_raactype:
                for kt in tmp_ktuple:
                    for gl in tmp_glambda:
                        yield self._encode(subtype=st, raactype=rt, ktuple=kt, glambda=gl)

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]], str, int, int, int], Tuple[List[str], List[float]]],
            aa_groups: Dict[int, List[str]],
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        self.func = func
        self.subtype = subtype
        self.raactype = raactype
        self.ktuple = ktuple
        self.glambda = glambda
        _, tmp_raactype, tmp_ktuple, _ = self._make_iterables()
        if not all([i in aa_groups for i in tmp_raactype]):
            raise ValueError(self.__class__.__name__ + ": raactype value is not correct.")
        if not all([i in [1, 2, 3] for i in tmp_ktuple]):
            raise ValueError(self.__class__.__name__ + ": invalid ktuple.")


class PseKRAACType1Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t1.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t1.type1, t1.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType2Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t2.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t2.type1, t2.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType3AEncoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t3A.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t3A.type1, t3A.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType3BEncoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t3B.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t3B.type1, t3B.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType4Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t4.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t4.type1, t4.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType5Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t5.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t5.type1, t5.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType6AEncoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t6A.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t6A.type1, t6A.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType6BEncoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t6B.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t6B.type1, t6B.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType6CEncoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t6C.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t6C.type1, t6C.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType7Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t7.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t7.type1, t7.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType8Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t8.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t8.type1, t8.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType9Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t9.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t9.type1, t9.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType10Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t10.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t10.type1, t10.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType11Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t11.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t11.type1, t11.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType12Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t12.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t12.type1, t12.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType13Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t13.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t13.type1, t13.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType14Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t14.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t14.type1, t14.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType15Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t15.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t15.type1, t15.AAGroup, subtype, raactype, ktuple, glambda)


class PseKRAACType16Encoder(_PseKRAACEncoder):

    @staticmethod
    def get_valid_aa_groups() -> Dict[int, List[str]]:
        return t16.AAGroup

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            subtype: Union[str, Iterable[str]],
            raactype: Union[int, Iterable[int]],
            ktuple: Union[int, Iterable[int]],
            glambda: Union[int, Iterable[int]]):
        super().__init__(in_data, cores, t16.type1, t16.AAGroup, subtype, raactype, ktuple, glambda)
