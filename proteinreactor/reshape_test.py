import json

seqs_str = \
""">Seq_1
RSAPAAAI
>Seq_2
SAPAAAIA
>Seq_3
APAAAIAA
>Seq_4
PAAAIAAR
>Seq_5
AAAIAARV
>Seq_6
AAIAARVA
>Seq_7
AIAARVAG
>Seq_8
IAARVAGQ
>Seq_9
AARVAGQT
>Seq_10
ARVAGQTR
>Seq_11
RVAGQTRN
>Seq_12
VAGQTRNI
>Seq_13
AGQTRNIT
>Seq_14
GQTRNITV
>Seq_15
QTRNITVD
>Seq_16
TRNITVDP
>Seq_17
RNITVDPR
>Seq_18
NITVDPRL
>Seq_19
ITVDPRLF
>Seq_20
TVDPRLFK
>Seq_21
VDPRLFKK"""

names, seqs = [], []

for l in seqs_str.split("\n"):
    l = l.rstrip()
    if l.startswith(">"):
        names += [l]
    else:
        seqs += [l]

import numpy as np

preds = np.random.binomial(1, 0.3, 21)

orig_seq = "".join([seqs[0]] + [i[-1] for i in seqs[1:]])
print(orig_seq)
print("".join([str(c) for c in preds]))

idxs_all = np.array(range(len(preds)))
print({
    "sequence": orig_seq,
    "difficult_couplings": list(map(lambda x: int(x), idxs_all[preds > 0])),
    "easy_couplings": list(map(lambda x: int(x), idxs_all[preds == 0]))
})