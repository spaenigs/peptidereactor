from modlamp.core import read_fasta
import pandas as pd
import os, sys, re
sys.path.append(os.getcwd())
from nodes.utils.tertiary_structure_prediction.scripts.utils import *

TOKEN = config["token"]

rule motif_search:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta"
    output:
         f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv"
    run:
         seqs, _ = read_fasta(str(input))
         query = seqs[0]
         query_len = len(query)

         df_res = pd.DataFrame()
         for pdb_id in pdb_motif_search(motif=query):
             seqs, names = pdb_get_fasta(pdb_id)
             chain_id = "A"
             for seq, name in zip(seqs, names):
                 if query in seq:
                     chain_id = re.search("^\w{4}:(\w)|.*", name).group(1)
                     seq_from_fasta = seq
                     break
             seq_from_fasta = get_seq_from_fasta(pdb_id, chain_id)
             idx_start = seq_from_fasta.find(query)
             best_hit = f"{pdb_id}_{chain_id}"
             df_tmp = pd.DataFrame(
                 {"qseqid": [query], "sacc": [best_hit], "sstart": [idx_start], "send": [idx_start + query_len], "evalue": [-1]})
             df_res = pd.concat([df_res, df_tmp])

         df_res.to_csv(str(output))