from glob import glob
from modlamp.core import read_fasta

import pandas as pd

import os

df = pd.DataFrame([(d.replace("seqs.fasta", ""), len(read_fasta(d)[1])) for d in glob("data/*/seqs.fasta")])
df.columns = ["path", "nr_of_seqs"]
df["present"] = df["path"].apply(lambda p: os.path.exists(p + "pdb/"))
df.sort_values(by=["nr_of_seqs", "present"], inplace=True)
df.reset_index(inplace=True)

print(df)

