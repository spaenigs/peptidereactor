import os
from io import StringIO

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from modlamp.core import read_fasta, save_fasta

import joblib as jl
import pandas as pd

TOKEN = config["token"]
TARGET_PDBS = config["pdbs_out"]
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
         temp(f"data/temp/{TOKEN}/{{seq_name}}.tmp.fasta")
    run:
         seq_tuples, seq_class = jl.load(str(input))
         save_fasta(str(output), sequences=[seq_tuples[0][1]], names=[seq_tuples[0][0]])

rule distribute:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.tmp.fasta",
         "peptidereactor/db/swiss_port/proteins.db"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta"
    run:
         seqs, names = read_fasta(input[0])
         seq, name = seqs[0], names[0]

         # TODO sort by subject total protein length

         if len(seq) < 30:
             header = ["qseqid", "sacc", "sstart", "send", "evalue", "qseq", "sseq"]
             cline = NcbiblastpCommandline(
                 task="blastp-short" if len(seq) < 15 else "blastp",
                 db=input[1],
                 query=input[0],
                 evalue=200000, # include distinct hits
                 max_target_seqs=1,
                 max_hsps=1,
                 outfmt="10 " + " ".join(header))
             stdout, stderr = cline()

             df_res = pd.read_csv(StringIO(stdout), names=header)

             if df_res.empty:
                 shell("touch {output[0]}")

             else:
                 id = df_res.loc[0, "sacc"]

                 fasta_lines = \
                     os.popen(f"blastdbcmd -db peptidereactor/db/pdb/in_structure/pdb.db -entry {id}").read()

                 if fasta_lines == "":
                     shell("touch {output[0]}")

                 else:
                     fasta = ""
                     for l in fasta_lines:
                         fasta += l
                     with open(output[0], "w") as f:
                         sequences = SeqIO.parse(StringIO(fasta), "fasta")
                         SeqIO.write(sequences, f, "fasta")

         else:
             save_fasta(output[0], seqs, names)

rule build_feature:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         expand("peptidereactor/RaptorX/databases/{database}/{database}.moved.txt",
                database=["nr70", "nr90"]),
         "peptidereactor/RaptorX/setup.pl"
    output:
         PROFILE_DIR + "{seq_name}.tgt"
    shell:
         """
         touch {output};
         export OLDWD=$PWD; cd peptidereactor/RaptorX/;
         ./buildFeature -i $OLDWD/{input[0]} -o $OLDWD/{output} -c 1 || \
            echo "No alignments found for {wildcards.seq_name}";
         cd - 1> /dev/null;
         """


