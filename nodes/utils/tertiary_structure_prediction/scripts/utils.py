import requests
from Bio import SeqIO, BiopythonWarning
from Bio.PDB import PDBList, PDBParser, MMCIFParser
from Bio.SeqUtils import seq1
import os
import re
from io import StringIO
import pandas as pd
from modlamp.core import read_fasta
import warnings


warnings.simplefilter('ignore', BiopythonWarning)


def get_structure(pdb_id):
    pdb_dir = "data/temp/pdbs/"
    # try to download as pdb file first
    pdb_path = PDBList()\
        .retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format="pdb")
    if os.path.exists(pdb_path):
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_path)
    else:
        # if fails, try cif format (default format)
        cif_path = PDBList() \
            .retrieve_pdb_file(pdb_id, pdir=pdb_dir)
        parser = MMCIFParser()
        structure = parser.get_structure(pdb_id, cif_path)
    return structure


def get_chain(structure, chain_id):
    return [chain for chain in list(structure.get_models())[0] if chain.get_id() == chain_id][0]


def get_seq_from_fasta(pdb_id, chain_id):
    print((pdb_id, chain_id))
    fasta_handle = os.popen(f"blastdbcmd -db peptidereactor/db/pdb/pdb.db -entry {pdb_id}") \
        .read().rstrip()
    return \
        [str(r.seq) for r in SeqIO.parse(StringIO(fasta_handle), "fasta")
         if f"_{chain_id}" in r.id][0]


def get_seq_from_pdb(chain):
    seq_from_pdb = seq1("".join([residue.get_resname() for residue in chain]))
    seq_from_pdb_ics = [residue.get_id()[1] for residue in chain]
    return seq_from_pdb, seq_from_pdb_ics


def pdb_motif_search(motif):
    query_text = f"""
<?xml version="1.0" encoding="UTF-8"?>
<orgPdbQuery>
<queryType>org.pdb.query.simple.MotifQuery</queryType>
<description>Motif Query For: {motif}</description>
<motif>{motif}</motif>
</orgPdbQuery>
"""
    # get best hits with descending order
    url = 'http://www.rcsb.org/pdb/rest/search/?sortfield=rank%20Descending'
    result = requests.post(
        url, data=query_text,
        headers={'Content-Type': 'application/x-www-form-urlencoded'}
    ).text
    ids = [t.split(":")[0] for t in result.rstrip().split("\n")]
    # in case of no hits
    ids = [] if ids[0] == "null" else ids
    return ids[:5] if len(ids) >= 5 else ids


def pdb_get_fasta(id):
    url = f"https://www.rcsb.org/pdb/download/downloadFastaFiles.do?" + \
          f"structureIdList={id}&compressionType=uncompressed"
    response = requests.get(url)
    fasta_handle = response.content.decode()
    seqs, names = [], []
    for r in SeqIO.parse(StringIO(fasta_handle), "fasta"):
        seqs += [str(r.seq)]
        names += [r.id]
    return seqs, names


def get_candidate_structures(path_to_fasta, path_to_hits):
    seqs, names = read_fasta(path_to_fasta)
    query = seqs[0]
    query_len = len(query)
    df = pd.read_csv(path_to_hits)
    print(df)
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
            {"pdb_id": [pdb_id], "chain_id": [chain_id], "best_hit": [best_hit], "start": [start], "end": [end]})
        df_res = pd.concat([df_res, df_tmp])
    return df_res
