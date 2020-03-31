from modlamp.core import read_fasta
import pandas as pd
import yaml
from more_itertools import flatten
import re

TOKEN = config["token"]

rule all:
    input:
         config["csv_out"]

rule encode:
    input:
         config["fasta_in"],
         config["profile"]
    output:
         temp(f"data/temp/{TOKEN}/enco_seqs.yaml")
    run:
         seqs, names = read_fasta(str(input[0]))
         fastas = [[n, s] for s, n in zip(seqs, names)]

         # adapted from peptidereactor/iFeature/codes/DisorderB.py
         enco = {"enco_seqs": {}}
         for name, sequence in fastas:

             with open(str(input[1]) + f"{name}.dis") as f:
                 records = f.readlines()

             tag = 0
             for i in range(len(records)):
                 if re.search('^-------', records[i]):
                     tag = i
                     break
             records = records[tag+1:-1]

             disValue = []
             for line in records:
                 array = line.rstrip().split() if line.rstrip() != '' else None
                 key = array[3] if array[3] == 'D' else 'O'
                 disValue.append(key)

             enco["enco_seqs"][name] = {"encoded_seq": disValue}

         with open(str(output), mode="w") as f:
            yaml.safe_dump(enco, f)

rule make_sliding_windows:
    input:
         config["fasta_in"],
         config["classes_idx_in"],
         config["classes_in"],
         f"data/temp/{TOKEN}/enco_seqs.yaml"
    output:
         config["csv_out"]
    run:
         seqs, names = read_fasta(str(input[0]))

         with open(str(input[1])) as f1, \
                 open(str(input[2])) as f2, \
                 open(str(input[3])) as f3:
             classes_idx = list(map(lambda l: int(l.rstrip()), f1.readlines()))
             windowed_classes = yaml.safe_load(f2)
             encoded_seqs = yaml.safe_load(f3)

         df_res, classes_res = pd.DataFrame(), []
         for seq, name, class_idx in zip(seqs, names, classes_idx):
             values = list(windowed_classes[class_idx].values())[0]
             for i, v in enumerate(values, start=1):
                 start, end = v["range"]
                 classes_res += [v["class"]]
                 encoded_seq = \
                     encoded_seqs["enco_seqs"][name]["encoded_seq"]
                 encoded_seq_window = \
                     [{'D':[0, 1], 'O':[1, 0]}[sse]
                      for sse in encoded_seq[start:end]]
                 df_tmp = pd.DataFrame({f"{name}_part_{str(i)}": flatten(encoded_seq_window)})\
                     .transpose()
                 df_res = pd.concat([df_res, df_tmp])

         df_res["y"] = classes_res
         df_res.to_csv(str(output))