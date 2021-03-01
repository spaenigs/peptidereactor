from modlamp.core import read_fasta
from scipy.stats import describe
from glob import glob

import numpy as np
import pandas as pd

lens_pos, lens_neg = [], []
data = []
for p in sorted(glob("data/*/seqs.fasta")):
    classes_path = p.replace("seqs.fasta", "classes.txt")
    with open(classes_path) as f:
        classes = [int(l.rstrip()) for l in f.readlines()]
    seqs, names = read_fasta(p)
    zipped = dict(zip(seqs, classes))
    len_pos = [len(seq) for seq, class_ in zipped.items() if class_ == 1]
    lens_pos += len_pos
    len_neg = [len(seq) for seq, class_ in zipped.items() if class_ == 0]
    lens_neg += len_neg
    des_pos, des_neg = describe(len_pos), describe(len_neg)
    data += [{
        "dataset": p.split("/")[-2].replace("_", "\\_"),
        "pos_neg_nobs": f"{des_pos.nobs}, {des_neg.nobs}",
        "pos_minmax": f"{des_pos.minmax[0]}, {des_pos.minmax[1]}",
        "neg_minmax": f"{des_neg.minmax[0]}, {des_neg.minmax[1]}",
        "pos_mean_std": f"{np.round(des_pos.mean, 2)}, {np.round(np.sqrt(des_pos.variance), 2)}",
        "neg_mean_std": f"{np.round(des_neg.mean, 2)}, {np.round(np.sqrt(des_neg.variance), 2)}",
        "pos_median": np.median(len_pos),
        "neg_median": np.median(len_neg)
    }]

des_pos = describe(lens_pos)
des_neg = describe(lens_neg)
des = describe(lens_pos + lens_neg)

print("+", des_pos, np.median(lens_pos), np.sqrt(des_pos.variance))
print("-", des_neg, np.median(lens_neg), np.sqrt(des_neg.variance))
print(" ", des, np.median(lens_pos + lens_neg), np.sqrt(des.variance))

pd.DataFrame(data).to_csv("data.csv", sep="&", index=False, line_terminator="\\\\ \\hline \n")
