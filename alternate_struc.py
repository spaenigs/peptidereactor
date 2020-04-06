from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.PDB import PDBParser, PDBList, Dice
from Bio import BiopythonWarning, SeqIO
from io import StringIO
from Bio.SeqUtils import seq1
from modlamp.core import read_fasta
import pandas as pd
import re
import os
import warnings

warnings.simplefilter('ignore', BiopythonWarning)


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


def myIter(rows, query):

    size = len(rows)

    if size > 5:
        rows = rows[:5]

    # TODO case (3), that is, no hits are found

    blast_query = query
    query_len = len(query)

    if size == 0:
        print(len(os.popen(f"cat pdb_seqres.txt | grep {query} -B 1 | grep '>'").read().rstrip().split("\n")))
        for line in os.popen(f"cat pdb_seqres.txt | grep {query} -B 1 | grep '>'").read().rstrip().split("\n"):
            print(line)
            pdb_id, chain_id = re.search(">(\w{4})_(\w).*", line).groups()
            structure = get_structure(pdb_id)
            chain = get_chain(structure, chain_id)
            seq_from_pdb, seq_from_pdb_ics = get_seq_from_pdb(chain)
            idx_start = seq_from_pdb.find(blast_query)
            best_hit = blast_query
            start = seq_from_pdb_ics[idx_start]
            end = start + query_len - 1
            # the desired section is missing in the pdb:
            if idx_start == -1:
                start, end = -1, -1
                best_hit = "X" * query_len
            yield pdb_id, chain_id, best_hit, start, end

    for name, series in rows:
        # TODO sequence part is missing
        pdb_id, chain_id = series['sacc'].split("_")
        blast_start, blast_end = series[["sstart", "send"]]

        structure = get_structure(pdb_id)
        chain = get_chain(structure, chain_id)

        seq_from_fasta = get_seq_from_fasta(pdb_id, chain_id)
        seq_from_pdb, seq_from_pdb_ics = get_seq_from_pdb(chain)

        blast_hit = seq_from_fasta[blast_start:blast_end]
        query_len = len(query)
        idx_start = seq_from_pdb.find(blast_query)
        best_hit = blast_query

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

        yield pdb_id, chain_id, best_hit, start, end


seqs, _ = read_fasta("short_query.fasta")
query = seqs[0]

header = ["qseqid", "sacc", "sstart", "send", "evalue"]

cline = NcbiblastpCommandline(
    task="blastp-short", db="db/pdb/pdb.db", query="short_query.fasta", outfmt="10 " + " ".join(header)
)

stdout, stderr = cline()

df = pd.read_csv(StringIO(stdout), names=header)

print(len(list(set(df["qseqid"]))))

for pdb_id, chain_id, hit, start, end in myIter(list(df.iterrows()), query):
    # TODO sort
    print((pdb_id, chain_id, hit, start, end))

# hits_found = df.shape[0]

# for name, series in df.iterrows():
# print("CCC")
# if hits_found == 0:
#     pdb_id, chain_id = "", ""
#     for line in os.popen("cat pdb_seqres.txt | grep AANDGPMP -A 2 | grep '>'").read().rstrip().split("\n"):
#         pdb_id, chain_id = re.search(">(\w{4})_(\w).*", line).groups()
#         break
# else:

# pdb_id, chain_id = series['sacc'].split("_")
#
# print(pdb_id)
#
# structure = get_structure(pdb_id)


# seq_pdb = [residue.get_resname() for residue in chain]
# seq_pdb_ics = [residue.get_id()[1] for residue in chain]
# seq_pdb_one_letter = seq1("".join(seq_pdb))
#
# idx_start = seq_pdb_one_letter.find(query)
#
# # if non-standard amino acids in sequence, use second-best blast hit
# if idx_start == -1 and hits_found > 1:
#     continue
# elif idx_start == -1 and hits_found == 1:
#     # find seqs with non-standard AAs
#     hits = []
#     for hit in re.finditer(f"[({query})X]{{{query_len}}}(?=[YEDINFGACMRHLSPWQKVT])", seq_pdb_one_letter):
#         hits += [dict(hit=hit.group(), idx_start=hit.start())]
#     # use section with least nr. of "X"
#     idx_start = sorted(hits, key=lambda d: len([c for c in d["hit"] if c == "X"]))[0]["idx_start"]
#
# print(f"Found query at {seq_pdb_ics[idx_start]} in chain '{chain_id}'")
#
# start = seq_pdb_ics[idx_start]
# end = start + query_len - 1

# start, end, is_valid = get_indices(structure, chain_id)
#
# if is_valid:
#     Dice.extract(structure, chain_id, start, end, f"{pdb_id}_{chain_id}_{start}-{end}.pdb")
#     break


# print((start,end))
#
#
#
# print(len([r for r in parser.header["missing_residues"] if r["chain"] == chain_id]))
#
# seqs = ""
# indices = []
# for i in [chain for chain in list(structure.get_models())[0] if chain.get_id() == "A"][0]:
#     res, idx = i.get_resname(), i.get_id()[1]
#     # print((res, idx))
#     seqs += res
#     indices += [idx]
#
#
#
# print("Query pos in pdb file: " + str(seq1(seqs).find("AANDGPMP")))
# print(seq3("AANDGPMP"))
# print(len(indices))
# print(parser.header)
