from abc import abstractmethod, ABC
from typing import List, Iterable, Union, Callable

from Bio.Align.Applications import MuscleCommandline

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import pandas as pd
import tempfile as tf
from subprocess import Popen, PIPE


class BaseEncoder(ABC):

    @staticmethod
    def get_meta_name(package: str, encoder: str, **kwargs):
        res = package + "_" + encoder
        for key, value in kwargs.items():
            res += "_{0}-{1}".format(key, value)
        return res

    @staticmethod
    def to_df(x: np.array, y: List[int], col_names: Iterable[str], row_names: List[str], meta_name: str) -> pd.DataFrame:
        df = pd.DataFrame(x)
        df.columns = col_names
        df.index = row_names
        df["y"] = y
        df.meta_name = meta_name
        return df

    @staticmethod
    def _sequences_sufficent_length(
            sequences: List[str],
            gap: Union[int, Iterable[int]],
            min_len: Callable[[int], int]) -> bool:
        def filter_len(seqs, _gap):
            return len(list(filter(lambda seq: len(seq) <= min_len(_gap), seqs)))
        if type(gap) == range:
            for g in gap:
                if filter_len(sequences, g) > 0:
                    return False
        else:
            if filter_len(sequences, gap) > 0:
                return False
        return True

    @staticmethod
    def run_muscle(in_data: (List[List[str]], List[int])) -> (int, List[List[str]]):
        seq_records = list(map(lambda t: SeqRecord(Seq(t[1]), description=t[0], id=t[0]), in_data[0]))
        tempfile_fasta_in = tf.NamedTemporaryFile(mode="w")
        tempfile_fasta_out = tf.NamedTemporaryFile(mode="w")
        SeqIO.write(seq_records, tempfile_fasta_in.name, "fasta")
        cline = MuscleCommandline(input=tempfile_fasta_in.name, out=tempfile_fasta_out.name)
        o, e = cline()
        fastas_aligned = [
            [str(record.name), str(record.seq)]
            for record in sorted(
                SeqIO.parse(tempfile_fasta_out.name, "fasta"),
                # key=lambda record: int(record.name.replace("Seq_", ""))
                key=lambda record: record.name
            )
        ]
        return len(fastas_aligned[0][1]), fastas_aligned

    def batch_interpolate(self, encoded_seqs: List[List[float]], interpolate_to: int) -> List[List[float]]:
        seqs_as_str = str(encoded_seqs)\
            .replace("[", "c(").replace("]", ")").replace("c", "list", 1)
        with tf.NamedTemporaryFile(mode="w") as f:
            f.write(seqs_as_str)
            f.flush()
            p = Popen(
                ['encoder/R_scripts/interpolate_seq.R', f.name, str(interpolate_to)],
                stdin=PIPE, stdout=PIPE, stderr=PIPE
            )
            output, _ = p.communicate()
            encoded_seqs = eval(output.decode())
        return encoded_seqs

    @abstractmethod
    def encode(self):
        pass

    def __init__(self, in_data: (List[List[str]], List[int]), cores: int):
        self.in_data = in_data
        self.cores = cores
        self.files = []


class BaseEncoderEqualLength(BaseEncoder):

    def _encode(self, func, **kwargs) -> Union[pd.DataFrame, Iterable[pd.DataFrame]]:
        resm = []
        for tup in self.in_data[0]:
            names, encoded_seq = func([tup])
            resm.append(encoded_seq[1:])
        interpolated_seqs = self.batch_interpolate(resm, self.interpolate_to)
        cn = map(
            lambda x: "{}.INTERP{}.{}".format(func.func.__name__, self.interpolate_to, x),
            range(self.interpolate_to)
        )
        return self.to_df(
            x=interpolated_seqs,
            y=self.in_data[1],
            col_names=cn,
            row_names=self.sequences,
            meta_name=self.get_meta_name(
                package="ifeature",
                encoder=type(self).__name__.lower(),
                **{**{"interpol": int(self.interpolate_to)}, **kwargs}
            )
        )

    def _encode_binary_with_gaps(self, fastas_aligned: List[List[str]], func, binary_len: int) -> pd.DataFrame:
        resm = []
        for tup in zip(fastas_aligned, self.in_data[0]):
            _, seq_aligned = tup[0]
            description, seq_original = tup[1]
            tripletts = np.array(
                func(fastas=[[description, seq_original]])[1][1:]
            ).reshape((-1, binary_len))
            widx = list(zip(seq_original, tripletts))
            tmp = list(reversed(widx.copy()))
            res = []
            for c in seq_aligned:
                if c == "-":
                    res.append([0] * binary_len)
                else:
                    triplett = tmp.pop()
                    res.append(list(triplett[1]))
            resm.append(list(np.array(res).ravel()))
        x = np.array(resm)
        cn = map(
            lambda _x: "{}.INTERP{}.{}".format(func.func.__name__, self.interpolate_to, _x),
            range(self.interpolate_to)
        )
        rn = list(map(lambda _tup: _tup[1], fastas_aligned))
        return self.to_df(
            x=x,
            y=self.in_data[1],
            col_names=cn,
            row_names=rn,
            meta_name=self.get_meta_name(package="ifeature", encoder=type(self).__name__.lower(), **{"interpol": int(self.interpolate_to / binary_len)})
        )

    def __init__(
            self,
            in_data: (List[List[str]], List[int]),
            cores: int,
            interpolation_factor: int):
        super().__init__(in_data, cores)
        self.sequences = [tup[1] for tup in self.in_data[0]]
        self.interpolate_to = int(np.median([len(seq) for seq in self.sequences])) * interpolation_factor
