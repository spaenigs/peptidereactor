from io import StringIO

import requests
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

import pandas as pd

df_res = pd.DataFrame()

for seq_len in range(50, 10000, 50)[:2]:

    print(seq_len)

    limit = 5000
    url = f"https://www.uniprot.org/uniprot/?query=length%3A[{seq_len}+TO+{seq_len}]&sort=score&format=fasta&limit={limit}"
    response = requests.get(url)
    fasta_handle = response.content.decode()

    with open("tmp.fasta", "w") as f:
        SeqIO.write(SeqIO.parse(StringIO(fasta_handle), "fasta"), f, "fasta")

    header = ["qseqid", "sacc", "sstart", "send", "evalue", "qseq", "sseq"]
    cline = NcbiblastpCommandline(
        db="peptidereactor/db/pdb/in_structure/pdb.db",
        query="tmp.fasta",
        num_threads=12,
        max_target_seqs=1,
        max_hsps=1,
        outfmt="10 " + " ".join(header))
    stdout, stderr = cline()

    df_blast = pd.read_csv(StringIO(stdout), names=header, index_col=0)

    df_res = pd.concat([df_res, df_blast])
    df_res.to_csv(f"blast_test_2.csv")

df_res.to_csv(f"blast_test_2.csv")
