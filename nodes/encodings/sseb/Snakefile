from modlamp.core import read_fasta
from iFeature import SSEB
import numpy as np
import pandas as pd

TOKEN = config["token"]

rule all:
    input:
         config["csv_out"]

rule encode:
    input:
         config["fasta_in"],
         config["profile"]
    output:
         f"data/temp/{TOKEN}/out.csv"
    run:
         def update_binary_encoded(_aligned_seq, _enco_seq):
             aligned_seq_len, binary_length = len(_aligned_seq), 3
             reshaped = np.array(_enco_seq).reshape((-1, binary_length))
             enco_seq_dict = dict(zip(range(aligned_seq_len), _aligned_seq))
             zeros_for_gaps = [[0]*binary_length] * aligned_seq_len
             res = dict(zip(range(aligned_seq_len), zeros_for_gaps))
             update_to = dict(zip([k for k, v in enco_seq_dict.items() if v != "-"], reshaped))
             res.update(update_to)
             return np.array(list(res.values()))\
                 .reshape((aligned_seq_len * binary_length, ))\
                 .tolist()

         seqs, names = read_fasta(str(input[0]))
         fastas = [[n, s] for s, n in zip(seqs, names)]

         enco = {}
         for tup in fastas:
             seq_name, aligned_seq = tup[0], tup[1]
             orig_seq = tup[1].replace("-", "")
             _, encoded_seq = SSEB.SSEB(fastas=[[seq_name, orig_seq]],
                                        path=str(input[1]))
             enco[seq_name] = update_binary_encoded(aligned_seq, encoded_seq[1:])

         pd.DataFrame(enco).transpose().to_csv(str(output))

rule dump:
    input:
         config["fasta_in"],
         config["classes_in"],
         f"data/temp/{TOKEN}/out.csv"
    output:
        config["csv_out"]
    run:
         seqs, names = read_fasta(str(input[0]))
         with open(str(input[1])) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         df = pd.read_csv(str(input[2]), index_col=0)
         df["y"] = -1

         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         for (name, (seq, class_)) in seq_tuples.items():
             df.loc[name, "y"] = class_

         df.to_csv(str(output))
