import pandas as pd
import yaml
from modlamp.core import read_fasta
from Bio.SeqUtils import seq3
from more_itertools import chunked

df_test = pd.read_csv("../data/pds/csv/non_empty_test/pds_window_length_11-waac_GEOR030106.csv", index_col=0)

with open("../data/temp/fba382c950e4/best_model_0_yprob.yaml") as f:
    probas = yaml.safe_load(f)

sequences, names = read_fasta("data/pds_window_length_11/seqs.fasta")
seq_tuples = dict((name, tup) for name, tup in zip(names, sequences))

with open("../data/pds/series.yaml") as f:
    series = yaml.safe_load(f)

df_test["probas"] = probas

df_test = df_test.loc[(df_test["probas"] >= 0.5) & (df_test["y"] == 0), :]

for seq_name in df_test.index:
    one_letter = seq_tuples[seq_name]
    three_letter = seq3(one_letter)
    three_letter = "-".join("".join(l) for l in chunked(three_letter, 3))
    for i in [f"{d['name']},{three_letter}" for d in series if three_letter in d["seq"]]:
        print(i)

# print(series)
