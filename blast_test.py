from io import StringIO

import requests
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

import pandas as pd

df_res = pd.DataFrame()


def download_proteins(fasta="random_proteins.fasta"):
    for seq_len in range(50, 10000, 50):
        print(seq_len)
        limit = 1000
        url = f"https://www.uniprot.org/uniprot/?query=length%3A[{seq_len}+TO+{seq_len}]&format=fasta&limit={limit}"
        response = requests.get(url)
        fasta_handle = response.content.decode()
        with open(fasta, "a") as f:
            SeqIO.write(SeqIO.parse(StringIO(fasta_handle), "fasta"), f, "fasta")


def run_blast(fasta="random_proteins.fasta"):
    header = ["qseqid", "sacc", "sstart", "send", "evalue", "qseq", "sseq", "qlen"]
    cline = NcbiblastpCommandline(
        db="peptidereactor/db/pdb/in_structure/pdb.db",
        query=fasta,
        num_threads=28,
        max_target_seqs=1,
        max_hsps=1,
        outfmt="10 " + " ".join(header))
    stdout, stderr = cline()
    df_res = pd.read_csv(StringIO(stdout), names=header, index_col=0)
    df_res["length_2"] = range(0, df_res.shape[0])
    df_res.to_csv(f"blast_test.csv")


download_proteins()
# run_blast()
