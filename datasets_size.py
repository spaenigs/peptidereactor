from glob import glob
from modlamp.core import read_fasta

import pandas as pd

import os

df = pd.DataFrame([(d.replace("seqs.fasta", ""), len(read_fasta(d)[1])) for d in glob("data/*/seqs.fasta")])
df.columns = ["path", "nr_of_seqs"]
df["enc_pres"] = df["path"].apply(lambda p: int(os.path.exists(p + "csv/structure_based/")))
df["ben_pres"] = df["path"].apply(lambda p: int(os.path.exists(p + "benchmark/dataset_correlation.csv")))
df.sort_values(by="nr_of_seqs", inplace=True)
df.reset_index(drop=True, inplace=True)

print(df)
