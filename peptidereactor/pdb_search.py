import requests
from Bio.SeqUtils import seq1
from Bio.PDB import *
from modlamp.core import read_fasta

url = 'http://www.rcsb.org/pdb/rest/search'

get_query = lambda motif: f"""
<?xml version="1.0" encoding="UTF-8"?>
<orgPdbQuery>
<queryType>org.pdb.query.simple.MotifQuery</queryType>
<description>Motif Query For: {motif}</description>
<motif>{motif}</motif>
</orgPdbQuery>
"""

seqs, names = read_fasta("short_query.fasta")

cnt = 0
for motif in seqs:
    print(motif)
    queryText = get_query(motif)
    result = requests.post(url, data=queryText, headers={'Content-Type': 'application/x-www-form-urlencoded'}).text
    print(result)
    ids = [t.split(":")[0] for t in result.rstrip().split("\n")]
    print(ids)
    if ids[0] == "null":
        cnt += 1

print(cnt/len(seqs))

# pdbl = PDBList()
# path = pdbl.retrieve_pdb_file(ids[2], file_format="pdb")
#
# p = PDBParser()
# structure = p.get_structure('X', path)
#
# class ChainSelect(Select):
#     def accept_chain(self, chain):
#         print(f"XXX {chain.get_id()}")
#         if chain.get_id()==self.c:
#             return 1
#         else:
#             return 0
#
#     def __init__(self, c):
#         self.c = c
#
# for model in structure:
#     print(model)
#     for chain in model:
#         print(chain)
#         seq = []
#         for residue in chain:
#             seq += [residue.get_resname()]
#                 # for atom in residue:
#                 #     print(atom)
#         final_seq = "".join([seq1(s) for s in seq])
#         if motif in final_seq:
#             io = PDBIO()
#             io.set_structure(structure)
#             print(chain.get_id())
#             io.save(f"out_{chain.get_id()}.pdb", ChainSelect(chain.get_id()))
#             print(final_seq)
