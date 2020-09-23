import subprocess

from sklearn.ensemble import RandomForestClassifier
from imblearn import combine as co
from imblearn import pipeline as pl

from Bio.Data.IUPACData import protein_letters_3to1, protein_letters_1to3
from modlamp.core import save_fasta
from sklearn.metrics import f1_score
from snakemake import shell

import numpy as np
import pandas as pd
import altair as alt

import yaml
import re
import more_itertools

with open('data/bce_bcm/bc_data/series.yaml') as f:
    series_result = yaml.safe_load(f)

df_res = pd.DataFrame()

for k, sr in enumerate(series_result):

    try:
        shell("rm -r data/temp/bachem/*")
    except Exception:
        pass

    final_seqs, classes = [], []

    window = 7

    for d in series_result:

        seq = re.match("H-(.*?)-NH2", d["seq"]).group(1).split("-")

        if sr == d:
            print([protein_letters_3to1[aa] for aa in seq])
            continue

        # pos_diff = list(map(lambda ac_pos: int(ac_pos.replace("Ac-", "")) - 1, d["rf"].keys()))
        pos_diff = [int(ac_pos.replace("Ac-", "")) - 2 for ac_pos, ac in d["rf"].items() if ac >= 1.5]
        # if no product, remove sequence rest
        if d["product"] == 0.0:
            smallest_ac_pos = sorted(pos_diff)[0]
            # keep 'product' by reassign Ac-(x-1) -> product. Will be remove later anyway.
            idx = (smallest_ac_pos - 1)
            tmps, tmpp = seq, pos_diff
            rest, seq = seq[:idx], seq[idx:]
            pos_diff = list(map(lambda i: i - len(rest), pos_diff))
        is_diff_coupl = np.zeros(len(seq)).astype(int)
        # set position where difficult coupling occurs
        is_diff_coupl[pos_diff] = 1
        # skip iteration if sequence length is smaller than window size,
        # - 1 due to subsequent removal of 'product'
        if len(seq) - 1 < int(window):
            continue
        # 1) remove Ac-1 (product) from sequence
        # 2) get all windows from N- to C-terminal
        # 3) reverse windows, such that we can read them in  intuitive oder (left to right)
        # 4) reverse position of difficult couplings
        for part_seq, couplings in \
                zip(reversed(list(more_itertools.windowed(seq[1:], int(window)))),
                    reversed(list(more_itertools.windowed(is_diff_coupl[1:], int(window))))):
            # 5) reverse partial sequences, because synthesis takes palce from C- to N-terminus
            final_seqs += ["".join(map(lambda aa3: protein_letters_3to1[aa3], reversed(part_seq)))]
            classes += [list(reversed(couplings))[-1]]

    seq_names = list(map(lambda i: "Seq_" + str(i), range(1, len(final_seqs) + 1)))

    save_fasta("data/temp/bachem/seqs.fasta", sequences=final_seqs, names=seq_names)
    with open("data/temp/bachem/classes.txt", mode="a") as f:
        for c in classes:
            f.write(f"{c}\n")
            f.flush()

    # 1) encode sequences
    cmd = "./peptidereactor/run_pipeline -s bcm_september_final.smk --cores 12"
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    process.communicate()

    # 2) train model
    df = pd.read_csv("data/temp/bachem/csv/cksaagp/cksaagp_gap_2.csv", index_col=0)
    X, y = df.iloc[:, :-1].values, df["y"]
    RANDOM_STATE = 42
    pipeline = pl.make_pipeline(
        co.SMOTEENN(random_state=RANDOM_STATE),
        RandomForestClassifier(random_state=RANDOM_STATE)
    )
    pipeline.fit(X, y)

    # 3) encode test peptide
    final_seqs, classes = [], []

    d = sr
    seq = re.match("H-(.*?)-NH2", d["seq"]).group(1).split("-")

    # pos_diff = list(map(lambda ac_pos: int(ac_pos.replace("Ac-", "")) - 1, d["rf"].keys()))
    pos_diff = [int(ac_pos.replace("Ac-", "")) - 2 for ac_pos, ac in d["rf"].items() if ac >= 1.5]
    # if no product, remove sequence rest
    if d["product"] == 0.0:
        smallest_ac_pos = sorted(pos_diff)[0]
        # keep 'product' by reassign Ac-(x-1) -> product. Will be remove later anyway.
        idx = (smallest_ac_pos - 1)
        tmps, tmpp = seq, pos_diff
        rest, seq = seq[:idx], seq[idx:]
        pos_diff = list(map(lambda i: i - len(rest), pos_diff))
    is_diff_coupl = np.zeros(len(seq)).astype(int)
    # set position where difficult coupling occurs
    is_diff_coupl[pos_diff] = 1
    # skip iteration if sequence length is smaller than window size,
    # - 1 due to subsequent removal of 'product'
    if len(seq) - 1 < int(window):
        continue
    # 1) remove Ac-1 (product) from sequence
    # 2) get all windows from N- to C-terminal
    # 3) reverse windows, such that we can read them in  intuitive oder (left to right)
    # 4) reverse position of difficult couplings
    for part_seq, couplings in \
            zip(reversed(list(more_itertools.windowed(seq[1:], int(window)))),
                reversed(list(more_itertools.windowed(is_diff_coupl[1:], int(window))))):
        # 5) reverse partial sequences, because synthesis takes palce from C- to N-terminus
        final_seqs += ["".join(map(lambda aa3: protein_letters_3to1[aa3], reversed(part_seq)))]
        classes += [list(reversed(couplings))[-1]]

    seq_names = list(map(lambda i: "Seq_" + str(i), range(1, len(final_seqs) + 1)))

    save_fasta("data/temp/bachem/seqs.fasta", sequences=final_seqs, names=seq_names)
    with open("data/temp/bachem/classes.txt", mode="w") as f:
        for c in classes:
            f.write(f"{c}\n")
            f.flush()

    cmd = "./peptidereactor/run_pipeline -s bcm_september_final.smk --cores 12"
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    process.communicate()

    # 4) predict test peptide
    df = pd.read_csv("data/temp/bachem/csv/cksaagp/cksaagp_gap_2.csv", index_col=0)
    X_test, y_test = df.iloc[:, :-1].values, df["y"]
    y_pred_bal = pipeline.predict(X_test)

    print(f1_score(y_test, y_pred_bal))

    # 5) decode test peptide
    decoded_seq = []
    for i, s in enumerate(final_seqs):
        if i == 0:
            decoded_seq += [protein_letters_1to3[aa] for aa in s]
        else:
            decoded_seq += [protein_letters_1to3[s[-1]]]
    fseq = [seq[0]] + list(reversed(decoded_seq))
    fcla = [0] + list(reversed(y_test.to_list())) + [0, 0, 0, 0, 0, 0]
    fcla_pred = [0] + list(reversed(y_pred_bal)) + [0, 0, 0, 0, 0, 0]

    # print("-".join(fseq))
    # print(fcla)
    # print(fcla_pred)

    # 6) save in df for later vis
    xvals = [f"Ac_{i+2}" for i in range(len(fseq))]
    seq_name =d["name"]
    source = pd.concat([
        pd.DataFrame({
            "x": xvals,
            "y": [f"{d['name']} (lab)" for i in range(len(fseq))],
            "seq": fseq, "coupl": [">=1.5" if c == 1 else "<1.5" for c in fcla],
            "f1": [f1_score(y_test, y_pred_bal) for i in range(len(fseq))],
            "seq_name": [seq_name for i in range(len(fseq))]
        }),
        pd.DataFrame({
            "x": xvals,
            "y": [f"{d['name']} (pred)" for i in range(len(fseq))],
            "seq": fseq, "coupl": [">=1.5" if c == 1 else "<1.5" for c in fcla_pred],
            "f1": [f1_score(y_test, y_pred_bal) for i in range(len(fseq))],
            "seq_name": [seq_name for i in range(len(fseq))]
        })
    ])

    df_res = pd.concat([df_res, source])

    # if k == 5:
    #     break

df_res.to_csv("source.csv")