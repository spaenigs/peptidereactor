from Bio.PDB import Dice
from modlamp.core import save_fasta

import os, sys
sys.path.append(os.getcwd())

from nodes.utils.tertiary_structure_prediction.scripts.utils import *

include:
    "setup_pdb.smk"

include:
    "download_pdb.smk"

TOKEN = config["token"]
TARGET_PDBS = config["pdbs_out"]
TARGET_DIR = os.path.dirname(TARGET_PDBS[0]) + "/" \
    if type(TARGET_PDBS) == list else os.path.dirname(TARGET_PDBS) + "/"

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

rule motif_search:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         "peptidereactor/db/cifs/"
    output:
         f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv"
    run:
         seqs, _ = read_fasta(str(input[0]))
         query = seqs[0]
         query_len = len(query)

         df_res = pd.DataFrame()
         for pdb_id in pdb_motif_search(motif=query):

             seqs, names = pdb_get_fasta(pdb_id)
             chain_id = "A"

             for seq, name in zip(seqs, names):
                 if query in seq:
                     chain_id = re.search("^\w{4}:(\w{1,2})|.*", name).group(1)
                     seq_from_fasta = seq
                     break

             seq_from_fasta = get_seq_from_fasta(pdb_id, chain_id)
             # +1: adapt to blast output
             idx_start = seq_from_fasta.find(query) + 1
             best_hit = f"{pdb_id}_{chain_id}"

             df_tmp = pd.DataFrame(
                 {"qseqid": [query], "sacc": [best_hit], "sstart": [idx_start],
                  "send": [idx_start + query_len], "evalue": [-1]})

             df_res = pd.concat([df_res, df_tmp])

         df_res.to_csv(str(output))

rule canditate_structures:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv"
    output:
         f"data/temp/{TOKEN}/structure_candidates_{{seq_name}}.csv"
    # threads:
    #     1000
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
         df = pd.read_csv(str(input), dtype={"pdb_id":str})

         if df.empty:
            shell("touch {output}")
         else:
             pdb_id, chain_id, start, end = \
                 df.loc[0, ["pdb_id", "chain_id", "start", "end"]]
             print((pdb_id, chain_id, start, end))
             structure = get_structure(pdb_id)
             Dice.extract(structure, chain_id, start, end, str(output))