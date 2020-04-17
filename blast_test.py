from io import StringIO

import requests
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

import pandas as pd

df_res = pd.DataFrame()
for seq_len in range(50, 10000, 50):
    print(seq_len)
    url = f"https://www.uniprot.org/uniprot/?query=length%3A[{seq_len}+TO+{seq_len}]&sort=score&format=fasta&limit=5000"
    response = requests.get(url)
    fasta_handle = response.content.decode()
    lens = []
    no_hits = 0
    for r in SeqIO.parse(StringIO(fasta_handle), "fasta"):
        with open("tmp.fasta", "w") as f:
            f.write(f">{r.id}\n{str(r.seq)}")
        header = ["evalue"] # ["qseqid", "sacc", "sstart", "send", "evalue", "qseq", "sseq"]
        cline = NcbiblastpCommandline(
            db="peptidereactor/db/pdb/in_structure/pdb.db",
            query="tmp.fasta",
            outfmt="10 " + " ".join(header))
        stdout, stderr = cline()
        df_blast = pd.read_csv(StringIO(stdout), names=header)
        if df_blast.shape[0] > 0:
            lens += [df_blast.evalue.values[0]]
        else:
            no_hits += 1
    df_tmp = pd.DataFrame({"evalues": lens, "no_hit": no_hits})
    df_tmp.index = [f"len_{seq_len}_{i}" for i  in df_tmp.index]
    df_res = pd.concat([df_res, df_tmp])
    df_res.to_csv("blast_test.csv")

df_res.to_csv("blast_test.csv")

