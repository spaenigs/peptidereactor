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
               seq_name=read_fasta(config["fasta_in"])[1][:20])

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
         f"data/temp/{TOKEN}/{{seq_name}}.fasta"
    run:
         seq_tuples, seq_class = jl.load(str(input))
         save_fasta(str(output), sequences=[seq_tuples[0][1]], names=[seq_tuples[0][0]])

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

# rule combine:
#     input:
#          f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv",
#          f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv"
#     output:
#          f"data/temp/{TOKEN}/combined_result_{{seq_name}}.csv"
#     run:
#          df_blast = pd.read_csv(str(input[0]))
#          df_motse = pd.read_csv(str(input[1]))
#
#          blast_hits = df_blast.shape[0]
#          motse_hits = df_motse.shape[0]
#
#          if blast_hits >= 5:
#              df_blast.iloc[:5, :].to_csv(str(output))
#          elif 0 < blast_hits < 5:
#              df_blast.to_csv(str(output))
#          elif motse_hits >= 5:
#              df_motse.iloc[:5, :].to_csv(str(output))
#          elif 0 < motse_hits < 5:
#              df_motse.iloc[:5, :].to_csv(str(output))
#          else:
#              pd.DataFrame({"qseqid": [], "sacc": [], "sstart": [], "send": [], "evalue": []})\
#                  .to_csv(str(output))

rule canditate_structures:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         # f"data/temp/{TOKEN}/combined_result_{{seq_name}}.csv"
         f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv"
    output:
         f"data/temp/{TOKEN}/structure_candidates_{{seq_name}}.csv"
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

# rule collect:
#     input:
#          expand(TARGET_DIR + "{seq_name}.pdb",
#                 seq_name=read_fasta(config["fasta_in"])[1])
#     output:
#          directory(TARGET_DIR)
#     run:
#          for path in list(input):
#              shell(f"cp {path} {output}")
