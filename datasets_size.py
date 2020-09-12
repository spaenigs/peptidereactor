from glob import glob
from modlamp.core import read_fasta

import pandas as pd

df = pd.DataFrame([(d, len(read_fasta(d)[1])) for d in glob("data/*/seqs.fasta")])
df.columns = ["path", "nr_of_seqs"]
df.sort_values(by="nr_of_seqs", inplace=True)
df.reset_index(inplace=True)

print(df)

