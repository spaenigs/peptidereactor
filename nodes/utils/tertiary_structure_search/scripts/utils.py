from Bio import BiopythonWarning, SeqIO
from Bio.PDB import MMCIFParser, Dice, PDBParser
from Bio.SeqUtils import seq1

import time
import requests
import re
import warnings

warnings.simplefilter('ignore', BiopythonWarning)


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


def get_seq_names(path_to_fasta):
    values = list(zip(*[(str(record.seq), record.id)
                        for record in SeqIO.parse(path_to_fasta, "fasta")]))
    if len(values) == 0:
        return []
    else:
        _, names = values
        return names


class Cif:

    def get_chain(self):
        return [chain for chain in list(self.structure.get_models())[0]
                if chain.get_id() == self.chain_id][0]

    def get_seq_from_pdb(self):
        seq_from_pdb = seq1("".join([residue.get_resname() for residue in self.chain]))
        seq_from_pdb = re.search("^X*(.*?)X*$", seq_from_pdb).group(1)
        seq_from_pdb_ics = [residue.get_id()[1] for residue in self.chain]
        return seq_from_pdb, seq_from_pdb_ics

    def dump_slice(self, motif, out_file):

        motif = motif.replace("-", "")
        start_on_indices = self.seq.find(motif)
        end_on_indices = start_on_indices + len(motif) - 1
        start, end = self.indices[start_on_indices], self.indices[end_on_indices]

        final_seq = \
            [r.get_resname() for r in self.chain.get_residues()
             if start <= r.get_id()[1] <= end]

        if "UNK" in final_seq:
            with open(out_file, "w") as f:
                f.write("")
                f.flush()
        else:
            Dice.extract(self.structure, self.chain_id, start, end, out_file)

    def __init__(self, pdb_id, chain_id, cif_dir, file_type="cif"):
        self.pdb_id = pdb_id
        self.chain_id = str(chain_id)
        if file_type == "cif":
            self.parser = MMCIFParser()
        else:
            self.parser = PDBParser()
        self.structure = self.parser.get_structure(pdb_id, cif_dir + f"{pdb_id}.{file_type}")
        self.chain = self.get_chain()
        self.seq, self.indices = self.get_seq_from_pdb()
