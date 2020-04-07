from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.PDB import Dice
from modlamp.core import save_fasta
import joblib as jl

import os, sys
sys.path.append(os.getcwd())
from nodes.utils.tertiary_structure_prediction.scripts.utils import *

TOKEN = config["token"]
TARGET_PDBS = config["pdbs_out"]
TARGET_DIR = os.path.dirname(TARGET_PDBS[0]) + "/" \
    if type(TARGET_PDBS) == list else os.path.dirname(TARGET_PDBS) + "/"

include:
    "setup_pdb.smk"

rule all:
    input:
        expand(TARGET_DIR + "{seq_name}.pdb",
               seq_name=read_fasta(config["fasta_in"])[1])

rule split_input_data:
    input:
         config["fasta_in"]
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta"
    run:
         seqs, names = read_fasta(str(input[0]))
         seq_tuples = dict((name, seq) for name, seq in zip(names, seqs))
         seq = seq_tuples[wildcards.seq_name]
         save_fasta(str(output), sequences=[seq], names=[wildcards.seq_name])

rule blast_search:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         "peptidereactor/db/pdb/pdb.db"
    output:
         f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv"
    priority:
        1000
    run:
         header = ["qseqid", "sacc", "sstart", "send", "evalue"]
         cline = NcbiblastpCommandline(
             task="blastp-short",
             db=str(input[1]),
             query=str(input[0]),
             outfmt="10 " + " ".join(header))
         stdout, stderr = cline()

         df_res = pd.read_csv(StringIO(stdout), names=header)
         blast_hits = df_res.shape[0]

         if blast_hits >= 5:
             df_res.iloc[:5, :].to_csv(str(output))
         else:
             df_res.to_csv(str(output))

rule canditate_structures:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv"
    output:
         f"data/temp/{TOKEN}/structure_candidates_{{seq_name}}.csv"
    threads:
        1000
    run:
         df_res = get_candidate_structures(
             path_to_fasta=str(input[0]), path_to_hits=str(input[1]))
         df_res.to_csv(str(output))

rule dump_cleaved_structure:
    input:
         f"data/temp/{TOKEN}/structure_candidates_{{seq_name}}.csv"
    output:
         TARGET_DIR + "{seq_name}.pdb"
    run:
         # TODO sort output
         df = pd.read_csv(str(input))

         if df.empty:
            shell("touch {output}")
         else:
             pdb_id, chain_id, start, end = \
                 df.loc[0, ["pdb_id", "chain_id", "start", "end"]]
             structure = get_structure(pdb_id)
             Dice.extract(structure, chain_id, start, end, str(output))