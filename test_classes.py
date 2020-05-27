from modlamp.core import read_fasta
from glob import glob

import pandas as pd

seqs, names = read_fasta("data/hiv_protease/seqs.fasta")
with open("data/hiv_protease/classes.txt") as f:
   classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))

for p in glob("data/hiv_protease/csv/structure_based/*.csv"):
   df = pd.read_csv(p, index_col=0)
   true, pred = [], []
   for n, s in df.iterrows():
      true += [s["y"]]
      pred += [seq_tuples[n][1]]

   if true != pred:
      print(p)


