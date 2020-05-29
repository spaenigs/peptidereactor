from io import StringIO
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from modlamp.core import read_fasta, save_fasta

import joblib as jl
import pandas as pd

import os

TOKEN = config["token"]
TARGET_DIR = config["target_dir"]
PROFILE_DIR = config["profile_dir"]

rule split_input_data:
    input:
         config["fasta_in"],
         config["classes_in"]
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.joblib")
    run:
         seqs, names = read_fasta(str(input[0]))
         with open(str(input[1])) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         seq_tuple = seq_tuples[wildcards.seq_name]
         jl.dump(value=([[wildcards.seq_name, seq_tuple[0]]], seq_tuple[1]), filename=str(output))

rule to_fasta:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.joblib"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.tmp")
    run:
         seq_tuples, seq_class = jl.load(str(input))
         save_fasta(str(output), sequences=[seq_tuples[0][1]], names=[seq_tuples[0][0]])

rule distribute:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.tmp",
         "peptidereactor/db/swiss_prot/proteins.db"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         f"data/temp/{TOKEN}/{{seq_name}}.csv"
    run:
         seqs, names = read_fasta(input[0])
         seq, name = seqs[0], names[0]

         if len(seq) < 30:
             header = ["qseqid", "sacc", "sstart", "send", "evalue", "qseq", "sseq"]
             cline = NcbiblastpCommandline(
                 task="blastp-short" if len(seq) < 15 else "blastp",
                 db=input[1],
                 query=input[0],
                 evalue=200000, # include distinct hits
                 outfmt="10 " + " ".join(header))
             stdout, stderr = cline()

             df_res = pd.read_csv(StringIO(stdout), names=header)

             if df_res.empty:
                 shell("touch {output[0]} {output[1]}")

             else:
                 df_res["length"] = -1
                 for idx in df_res.index:
                     id = df_res.loc[idx, "sacc"]
                     fasta_lines_tmp = \
                         os.popen(f"blastdbcmd -db {input[1]} -entry {id}").read()
                     for r in SeqIO.parse(StringIO(fasta_lines_tmp), "fasta"):
                         df_res.loc[idx, "length"] = len(r.seq)

                 df_res.sort_values(by="length", ascending=True, inplace=True)
                 df_res.to_csv(output[1])

                 id = df_res["sacc"].values[0]
                 fasta_lines = \
                         os.popen(f"blastdbcmd -db {input[1]} -entry {id}").read()

                 if fasta_lines == "":
                     shell("touch {output[0]}")

                 else:
                     seqs, names = [], []
                     for r in SeqIO.parse(StringIO(fasta_lines), "fasta"):
                         seqs += [str(r.seq)]
                         names += [f"{wildcards.seq_name}_{id}"]

                     save_fasta(output[0], seqs, names)

         else:
             save_fasta(output[0], seqs, names)
             pd.DataFrame().to_csv(output[1])

rule build_feature:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         expand("peptidereactor/RaptorX/databases/{database}/{database}.moved.txt",
                database=["nr70", "nr90"]),
         "peptidereactor/RaptorX/setup.pl"
    output:
         PROFILE_DIR + "{seq_name}.tgt",
         "peptidereactor/RaptorX/tmp/{seq_name}.diso",
         "peptidereactor/RaptorX/tmp/{seq_name}.psp"
    shell:
         """
         touch {output};
         export OLDWD=$PWD; cd peptidereactor/RaptorX/;
         ./buildFeature -i $OLDWD/{input[0]} -o $OLDWD/{output} -c 1 || \
            echo "No alignments found for {wildcards.seq_name}";
         cd - 1> /dev/null;
         """
