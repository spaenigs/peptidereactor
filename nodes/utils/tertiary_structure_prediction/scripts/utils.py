from Bio.PDB import MMCIFParser, Dice
from Bio.SeqUtils import seq1

import time
import requests
import re


def get_response(url):
    response = requests.get(url)
    cnt = 20
    while cnt != 0:
        if response.status_code == 200:
            return response.content.decode()
        else:
            time.sleep(1)
            cnt -= 1
    raise IOError(f"Some issues with PDB now. Try again later...\n(URL: {url}")


def get_chain(structure, chain_id):
    return [chain for chain in list(structure.get_models())[0] if chain.get_id() == chain_id][0]


def get_seq_from_pdb(chain):
    seq_from_pdb = seq1("".join([residue.get_resname() for residue in chain]))
    seq_from_pdb = re.search("^X*(.*?)X*$", seq_from_pdb).group(1)
    seq_from_pdb_ics = [residue.get_id()[1] for residue in chain]
    return seq_from_pdb, seq_from_pdb_ics


def dump_structure_slice(pdb_id, chain_id, motif, cif_dir, out_file):
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id, cif_dir + f"{pdb_id}.cif")
    chain = get_chain(structure, chain_id)
    seq, indices = get_seq_from_pdb(chain)
    start_on_indices = seq.find(motif)
    end_on_indices = start_on_indices + len(motif) - 1
    start, end = indices[start_on_indices], indices[end_on_indices]
    Dice.extract(structure, chain_id, start, end, out_file)