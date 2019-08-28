from typing import List, Callable, Tuple, Optional, Dict, Any, Iterable

import functools as ft
import numpy as np
import pandas as pd

from iFeature.codes.BINARY import BINARY
from iFeature.codes.BLOSUM62 import BLOSUM62
from iFeature.codes.ZSCALE import ZSCALE
from iFeature.codes.AAINDEX import AAINDEX

from iFeature.codes.AAC import AAC
from iFeature.codes.TPC import TPC
from iFeature.codes.DPC import DPC
from iFeature.codes.DDE import DDE
from iFeature.codes.GAAC import GAAC
from iFeature.codes.GDPC import GDPC
from iFeature.codes.GTPC import GTPC
from iFeature.codes.CTDC import CTDC
from iFeature.codes.CTDT import CTDT
from iFeature.codes.CTDD import CTDD

from iFeature.codes.checkFasta import checkFasta

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

from encoder.encoder import BaseEncoderEqualLength, BaseEncoder


class _ParamFreeEqualLengthEncoder(BaseEncoderEqualLength):

    def encode(self) -> pd.DataFrame:
        return self._encode(ft.partial(self.func))

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]]], Tuple[List[str], List[float]]],
            interpolation_factor: int):
        super().__init__(in_data, cores, interpolation_factor)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        self.func = func


class _ParamFreeDiffLengthEncoder(BaseEncoder):

    def _encode(self) -> pd.DataFrame:
        res = self.func(self.in_data[0], **dict({"order": None}))
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name("ifeature", type(self).__name__.lower())
        )

    def encode(self)-> pd.DataFrame:
        return self._encode()

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            func: Callable[[List[List[str]], Optional[Dict[str, Any]]], Tuple[List[str], List[float]]]):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        self.func = func


class Blosum62Encoder(_ParamFreeEqualLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, BLOSUM62, 20)


class ZscaleEncoder(_ParamFreeEqualLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, ZSCALE, 5)


class BinaryEncoder(_ParamFreeEqualLengthEncoder):

    # adapted from ifrom iFeature.codes.BINARY
    def _MY_BINARY(self, fastas, **kw):
        if not checkFasta(fastas):
            print('Error: for "BINARY" encoding, the input fasta sequences should be with equal length. \n\n')
            return 0

        AA = 'ARNDCQEGHILKMFPSTWYV-'
        encodings = []
        header = ['#']
        for i in range(1, len(fastas[0][1]) * 21 + 1):
            header.append('BINARY.F' + str(i))
        encodings.append(header)

        for i in fastas:
            name, sequence = i[0], i[1]
            code = [name]
            for aa in sequence:
                if aa == '-':
                    code = code + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
                    continue
                for aa1 in AA:
                    tag = 1 if aa == aa1 else 0
                    code.append(tag)
            encodings.append(code)
        return encodings

    def encode(self) -> pd.DataFrame:
        if self.run_msa:
            _, fastas_aligned = self.run_muscle(self.in_data)
            res = self._MY_BINARY(fastas_aligned)
        else:
            res = self._MY_BINARY(self.in_data[0])
        x = np.array(res[1:])[:, 1:]
        cn = res[0][1:]
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name("ifeature", type(self).__name__.lower())
        )

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, run_msa: bool):
        super().__init__(in_data, cores, BINARY, 1)
        self.run_msa = run_msa


class AAIndexEncoder(_ParamFreeEqualLengthEncoder):

    _aaindex_names = [
        'RADA880102', 'CIDH920104', 'NAKH920106', 'QIAN880131', 'NAKH920105', 'TANS770110', 'MCMT640101',
        'PONP800105', 'OOBM770105', 'MAXF760106', 'FAUJ880101', 'TANS770106', 'PALJ810113', 'ISOY800105',
        'JANJ790102', 'RICJ880112', 'GOLD730102', 'MEEJ800101', 'PALJ810108', 'CORJ870106', 'NADH010101',
        'PALJ810107', 'NADH010106', 'MUNV940101', 'GEOR030101', 'FUKS010110', 'CORJ870104', 'KANM800102',
        'LIFS790103', 'TANS770105', 'TAKK010101', 'RACS820111', 'TANS770107', 'DAYM780201', 'BASU050102',
        'QIAN880138', 'ZIMJ680104', 'BEGF750103', 'MONM990201', 'WILM950104', 'TSAJ990102', 'BEGF750101',
        'GEIM800111', 'ROBB760110', 'KUMS000102', 'HOPT810101', 'AURR980106', 'TANS770108', 'CEDJ970103',
        'CHOP780209', 'KLEP840101', 'OOBM770102', 'HOPA770101', 'BHAR880101', 'BUNA790103', 'JURD980101',
        'PALJ810106', 'ROBB760104', 'VASM830103', 'RICJ880106', 'KRIW790102', 'PONP800103', 'WEBA780101',
        'CHAM830101', 'GUOD860101', 'NADH010104', 'FASG890101', 'CEDJ970102', 'AURR980112', 'MIYS990103',
        'NISK800101', 'YUTK870101', 'OOBM770104', 'TANS770109', 'FAUJ880102', 'PONP800101', 'PALJ810105',
        'QIAN880120', 'MUNV940102', 'BLAM930101', 'EISD860103', 'RICJ880113', 'QIAN880103', 'WILM950102',
        'BROC820102', 'COSI940101', 'BASU050101', 'BIOV880102', 'WOLS870102', 'KUHL950101', 'CRAJ730101',
        'FODM020101', 'ISOY800107', 'ZHOH040102', 'QIAN880108', 'GEOR030107', 'NAKH900103', 'CHAM830107',
        'NADH010102', 'MIYS990101', 'DAWD720101', 'ENGD860101', 'SNEP660104', 'GEOR030105', 'PALJ810109',
        'ROSG850101', 'CHOP780211', 'GEOR030103', 'RADA880104', 'GRAR740102', 'QIAN880126', 'FAUJ880103',
        'RACS770101', 'FAUJ880108', 'BUNA790101', 'FUKS010106', 'GUYH850102', 'PRAM820102', 'MONM990101',
        'FUKS010101', 'VHEG790101', 'ONEK900102', 'CHAM830106', 'LEVM760106', 'FAUJ880107', 'GUYH850105',
        'NAKH900113', 'ARGP820103', 'LEVM780103', 'MAXF760102', 'LEVM760102', 'HUTJ700102', 'LEVM760104',
        'CHAM830105', 'ZIMJ680102', 'FUKS010108', 'JOND750101', 'LEVM760103', 'RICJ880117', 'SUYM030101',
        'GUYH850101', 'FUKS010107', 'RICJ880109', 'CHAM830102', 'CHOP780101', 'EISD860102', 'RACS820113',
        'LEVM780105', 'BIGC670101', 'QIAN880109', 'BURA740102', 'CHOC760102', 'RADA880105', 'PALJ810104',
        'ROBB760111', 'JANJ780102', 'FINA910104', 'CEDJ970104', 'TANS770101', 'ROSG850102', 'CHOC760104',
        'NAKH920108', 'TANS770104', 'CIDH920103', 'FINA770101', 'RICJ880102', 'RACS820105', 'PALJ810111',
        'CHOP780212', 'BEGF750102', 'KRIW790103', 'ROBB760106', 'KOEP990102', 'RACS820101', 'WOLS870101',
        'CHOP780205', 'MIYS850101', 'CHOP780208', 'YUTK870102', 'FAUJ880112', 'ISOY800103', 'BURA740101',
        'CEDJ970101', 'FAUJ880113', 'PTIO830101', 'VENT840101', 'NADH010107', 'QIAN880122', 'FASG760104',
        'PONP800102', 'QIAN880106', 'OOBM770103', 'WOEC730101', 'FAUJ880104', 'HUTJ700101', 'CHAM820101',
        'RICJ880104', 'AURR980110', 'MANP780101', 'MEEJ810102', 'NAKH920104', 'OOBM850104', 'FUKS010111',
        'BUNA790102', 'ZASB820101', 'QIAN880118', 'CHOP780215', 'EISD840101', 'PARJ860101', 'MIYS990102',
        'JANJ780101', 'RICJ880107', 'CHOP780204', 'AURR980109', 'ROSM880103', 'AURR980115', 'GEOR030104',
        'PALJ810101', 'GARJ730101', 'SUEM840102', 'OOBM850103', 'GEIM800110', 'CORJ870101', 'CIDH920102',
        'QIAN880115', 'WOLR790101', 'QIAN880110', 'CIDH920105', 'CHAM820102', 'PARS000102', 'MEEJ800102',
        'ROBB760112', 'VELV850101', 'QIAN880114', 'CASG920101', 'VINM940102', 'MEEJ810101', 'GRAR740103',
        'KUMS000101', 'AURR980107', 'RICJ880111', 'DESM900101', 'LEVM760101', 'TSAJ990101', 'OLSK800101',
        'BROC820101', 'CORJ870102', 'RACS820102', 'ISOY800102', 'KANM800103', 'NADH010105', 'KARP850101',
        'RICJ880108', 'PONJ960101', 'AURR980111', 'MEIH800102', 'JOND920102', 'PALJ810110', 'SNEP660101',
        'WOLR810101', 'QIAN880116', 'QIAN880124', 'MUNV940103', 'FAUJ880106', 'NAKH900109', 'KRIW790101',
        'SNEP660102', 'GOLD730101', 'FAUJ880105', 'NAKH900108', 'QIAN880132', 'LEVM780104', 'QIAN880119',
        'NAKH900110', 'FUKS010112', 'FAUJ880111', 'FAUJ880110', 'MAXF760101', 'WERD780104', 'CHOP780207',
        'GEOR030108', 'QIAN880136', 'AURR980104', 'JANJ780103', 'WOLS870103', 'TANS770102', 'AURR980103',
        'AURR980120', 'CHOC760103', 'FASG760102', 'FINA910103', 'ROBB790101', 'QIAN880104', 'JANJ790101',
        'PARS000101', 'CHAM830103', 'CORJ870107', 'CHOP780201', 'RACS820114', 'BULH740101', 'NOZY710101',
        'KARP850102', 'PTIO830102', 'QIAN880139', 'LIFS790101', 'RICJ880114', 'OOBM850105', 'LEVM780101',
        'PONP800104', 'GEIM800108', 'RICJ880103', 'BIOV880101', 'GEIM800106', 'ROBB760101', 'KANM800101',
        'FINA910101', 'GRAR740101', 'ISOY800106', 'LEVM780106', 'CHOP780214', 'GEIM800107', 'PUNT030101',
        'GEIM800102', 'GEIM800104', 'MIYS990105', 'GEIM800103', 'QIAN880123', 'BASU050103', 'NAGK730101',
        'QIAN880101', 'NAGK730102', 'ZHOH040101', 'COHE430101', 'FASG760105', 'YUTK870103', 'SUEM840101',
        'ROSM880101', 'RADA880108', 'JOND920101', 'QIAN880102', 'NAKH920107', 'JOND750102', 'BULH740102',
        'PALJ810114', 'SNEP660103', 'KRIW710101', 'AURR980116', 'RACS770102', 'WILM950103', 'QIAN880128',
        'TANS770103', 'ROBB760102', 'PLIV810101', 'HARY940101', 'MAXF760104', 'WIMW960101', 'GEIM800105',
        'CHOC760101', 'AURR980102', 'PONP800107', 'NAKH920103', 'QIAN880125', 'KOEP990101', 'EISD860101',
        'CHOP780202', 'SIMZ760101', 'OOBM770101', 'QIAN880112', 'RICJ880105', 'FUKS010109', 'AURR980105',
        'RADA880107', 'RACS820110', 'FUKS010104', 'ROBB760107', 'CHOP780216', 'FASG760101', 'AURR980119',
        'NAKH900102', 'ZIMJ680105', 'CHOC750101', 'COWR900101', 'QIAN880133', 'KHAG800101', 'RACS770103',
        'KANM800104', 'CHAM830108', 'PRAM900104', 'RACS820109', 'RACS820112', 'RICJ880101', 'PUNT030102',
        'NAKH920102', 'RACS820103', 'NAKH900112', 'ZIMJ680101', 'MUNV940105', 'ISOY800104', 'ROBB760108',
        'NAKH900104', 'DIGM050101', 'CORJ870105', 'PONP800106', 'KUMS000103', 'AURR980114', 'CEDJ970105',
        'KIDA850101', 'GEIM800101', 'VASM830102', 'MUNV940104', 'NAKH900111', 'FASG760103', 'ROSM880102',
        'NAKH920101', 'PALJ810103', 'AURR980113', 'ROBB760109', 'MEIH800101', 'CHOP780213', 'LEVM780102',
        'RADA880106', 'NAKH900107', 'GEOR030109', 'LEVM760107', 'ONEK900101', 'DESM900102', 'PALJ810115',
        'LEVM760105', 'GEOR030102', 'NAKH900105', 'VINM940101', 'WERD780102', 'VINM940104', 'JUNJ780101',
        'ROBB760113', 'LEWP710101', 'RADA880103', 'NAGK730103', 'RACS820108', 'ISOY800108', 'QIAN880111',
        'QIAN880127', 'PRAM900102', 'FUKS010102', 'CHOP780203', 'GUYH850104', 'CORJ870103', 'MIYS990104',
        'QIAN880113', 'RACS820106', 'ZHOH040103', 'RICJ880110', 'AURR980101', 'RICJ880116', 'HUTJ700103',
        'CHAM810101', 'FAUJ830101', 'PRAM820103', 'PALJ810112', 'GEOR030106', 'NAKH900106', 'FAUJ880109',
        'QIAN880130', 'ANDN920101', 'QIAN880105', 'VASM830101', 'FUKS010103', 'NISK860101', 'FUKS010105',
        'RICJ880115', 'CRAJ730102', 'MEIH800103', 'PONP930101', 'QIAN880134', 'FINA910102', 'ISOY800101',
        'AURR980108', 'MAXF760105', 'GEIM800109', 'JACR890101', 'QIAN880121', 'QIAN880137', 'CRAJ730103',
        'WERD780101', 'WILM950101', 'NADH010103', 'WERD780103', 'RACS820104', 'PRAM900101', 'LIFS790102',
        'SWER830101', 'NAKH900101', 'ROBB760105', 'PALJ810102', 'ROBB760103', 'ARGP820102', 'MAXF760103',
        'YUTK870104', 'QIAN880129', 'OOBM850101', 'AURR980117', 'VINM940103', 'CIDH920101', 'PRAM900103',
        'WARP780101', 'BAEK050101', 'JUKT750101', 'PRAM820101', 'KYTJ820101', 'OOBM850102', 'KARP850103',
        'PONP800108', 'CORJ870108', 'RADA880101', 'ARGP820101', 'PALJ810116', 'ZIMJ680103', 'AURR980118',
        'LAWE840101', 'RACS820107', 'KUMS000104', 'MITS020101', 'DAYM780101', 'CHAM830104', 'CHOP780210',
        'BLAS910101', 'QIAN880107', 'KIMC930101', 'QIAN880135', 'CHOP780206', 'QIAN880117'
    ]

    def encode(self) -> Iterable[pd.DataFrame]:
        for aaindex_name in self._aaindex_names:
            resm = []
            for tup in self.in_data[0]:
                names, vals = AAINDEX([tup])
                _, encoded_seq = list(
                    zip(
                        *filter(
                            lambda res_tup: aaindex_name in res_tup[0],
                            list(zip(names, vals))
                        )
                    )
                )
                resm.append(list(encoded_seq))
            interpolated_seqs = self.batch_interpolate(resm, self.interpolate_to)
            cn = map(
                lambda x: "{}.INTERP{}.{}".format(AAINDEX.__name__, self.interpolate_to, x),
                range(self.interpolate_to)
            )
            yield self.to_df(
                x=interpolated_seqs,
                y=self.in_data[1],
                col_names=cn,
                row_names=self.sequences,
                meta_name=self.get_meta_name("ifeature", type(self).__name__.lower(), **{"aaindex": aaindex_name, "interpol": self.interpolate_to})
            )

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int, index: str = None):
        super().__init__(in_data, cores, AAINDEX, 1)
        self._aaindex_names = self._aaindex_names if index is None else [index]


class AACEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, AAC)


class GAACEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, GAAC)


class GDPCEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, GDPC)


class TPCEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, TPC)


class DPCEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, DPC)


class DDEEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, DDE)


class GTPCEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, GTPC)


class CTDCEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, CTDC)


class CTDTEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, CTDT)


class CTDDEncoder(_ParamFreeDiffLengthEncoder):

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        super().__init__(in_data, cores, CTDD)
