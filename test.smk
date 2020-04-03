import os
import re
from io import StringIO

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils import seq1
from modlamp.core import read_fasta


def get_structure(pdb_id):
    # TODO handle case: Desired structure doesn't exists
    PDBList().retrieve_pdb_file(pdb_id, pdir="data/temp/pdbs/", file_format="pdb")
    parser = PDBParser()
    return parser.get_structure('PHA-L', f"data/temp/pdbs/pdb{pdb_id}.ent")


def get_chain(structure, chain_id):
    return [chain for chain in list(structure.get_models())[0] if chain.get_id() == chain_id][0]


def get_seq_from_fasta(pdb_id, chain_id):
    fasta_handle = os.popen(f"blastdbcmd -db db/pdb/pdb.db -entry {pdb_id}") \
        .read().rstrip()
    return \
        [str(r.seq) for r in SeqIO.parse(StringIO(fasta_handle), "fasta")
         if f"_{chain_id}" in r.id][0]


def get_seq_from_pdb(chain):
    seq_from_pdb = seq1("".join([residue.get_resname() for residue in chain]))
    seq_from_pdb_ics = [residue.get_id()[1] for residue in chain]
    return seq_from_pdb, seq_from_pdb_ics


rule blast:
    input:
         "short_query.fasta"
    output:
         "blast_result.csv"
    run:
         header = ["qseqid", "sacc", "sstart", "send", "evalue"]
         cline = NcbiblastpCommandline(
             task="blastp-short",
             db="db/pdb/pdb.db",
             query="short_query.fa",
             outfmt="10 " + " ".join(header))
         stdout, stderr = cline()
         pd.read_csv(StringIO(stdout), names=header).to_csv(str(output))

rule motif_search:
    input:
         "short_query.fasta"
    output:
         "motse_result.csv"
    run:
         seqs, names = read_fasta("short_query.fa")
         query = seqs[0]
         query_len = len(query)
         df_res = pd.read_csv()
         for line in os\
                 .popen(f"cat pdb_seqres.txt | grep {query} -B 1 | grep '>'").read().rstrip().split("\n"):
             pdb_id, chain_id = re.search(">(\w{4})_(\w).*", line).groups()
             seq_from_fasta = get_seq_from_fasta(pdb_id, chain_id)
             idx_start = seq_from_fasta.find(query)
             best_hit = query
             df_tmp = pd.DataFrame(
                 {"qseqid": query, "sacc": best_hit, "sstart": idx_start, "send": idx_start + query_len, "evalue": -1})
             df_res = pd.concat([df_res, df_tmp])
         df_res.to_csv(str(output))

rule combine:
    input:
         "blast_result.csv",
         "motse_result.csv"
    output:
         "combined_result.csv"
    run:
         df_blast = pd.read_csv(str(input[0]))
         df_motse = pd.read_csv(str(input[1]))

         blast_hits = df_blast.shape[0]
         motse_hits = df_motse.shape[0]

         if blast_hits >= 5:
             df_blast.iloc[:5, :].to_csv(str(output))
         elif 0 < blast_hits < 5:
             df_blast.to_csv(str(output))
         elif motse_hits >= 5:
             df_motse.iloc[:5, :].to_csv(str(output))
         elif 0 < motse_hits < 5:
             df_motse.iloc[:5, :].to_csv(str(output))
         else:
             pd.DataFrame().to_csv()

rule extract_indices:
    input:
         "short_query.fasta",
         "combined_result.csv"
    output:
         "structures.csv"
    run:
         seqs, names = read_fasta(str(input[0]))
         query = seqs[0]
         query_len = len(query)

         df = pd.read_csv(str(input[1]))

         df_res = pd.DataFrame()
         for name, series in df.iterrows():

             # TODO sequence part is missing
             pdb_id, chain_id = series['sacc'].split("_")
             blast_start, blast_end = series[["sstart", "send"]]

             structure = get_structure(pdb_id)
             chain = get_chain(structure, chain_id)

             seq_from_fasta = get_seq_from_fasta(pdb_id, chain_id)
             seq_from_pdb, seq_from_pdb_ics = get_seq_from_pdb(chain)

             blast_hit = seq_from_fasta[blast_start:blast_end]
             query_len = len(query)
             idx_start = seq_from_pdb.find(query)
             best_hit = query

             # if non-standard amino acids in sequence
             if idx_start == -1:
                 hits = []
                 for hit in re.finditer(f"[({query})X]{{{query_len}}}(?=[YEDINFGACMRHLSPWQKVT])",
                                        seq_from_pdb):
                     hits += [dict(hit=hit.group(), idx_start=hit.start())]
                 # use section with least nr. of "X"
                 if len(hits) > 0:
                     best_hit_dict = sorted(hits, key=lambda d: len([c for c in d["hit"] if c == "X"]))[0]
                     best_hit = best_hit_dict["hit"]
                     idx_start = best_hit_dict["idx_start"]

             start = seq_from_pdb_ics[idx_start]
             end = start + query_len - 1

             # the desired section is missing in the pdb:
             if idx_start == -1:
                 start, end = -1, -1
                 best_hit = "X" * query_len

             print((query, blast_hit, best_hit))

             df_tmp = pd.DataFrame(
                 {"pdb_id": pdb_id, "chain_id": chain_id, "best_hit": best_hit, "start": start, "end": end})
             df_res = pd.concat([df_res, df_tmp])

         df_res.to_csv(str(output))

rule get_structure:
    input:
         "structures.csv"
    output:
         "struc.pdb"
    run:
         # TODO sort output
         pass