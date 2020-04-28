import warnings
from io import StringIO
from glob import glob
from Bio import SeqIO, BiopythonWarning
from Bio.PDB import MMCIFParser, PDBList
from Bio.SeqUtils import seq1

import re

from nodes.utils.tertiary_structure_search.scripts.fetch_missing_ids \
    import fetch_and_save

warnings.simplefilter('ignore', BiopythonWarning)

rule all:
    input:
         "peptidereactor/db/pdb/in_structure/pdb.db"

# rule download_pdb:
#     output:
#          directory("peptidereactor/db/cifs/")
#     run:
#          attention: runs 3 to 4 days for whole pdb!
         # pl = PDBList(pdb=output[0])
         # pl.flat_tree = True

rule parse_cif_file:
    input:
         "peptidereactor/db/cifs/{id}.cif"
    output:
         "data/temp/fastas_from_cifs/{id}.fasta"
    run:
         id = wildcards.id

         attempt = 1
         while attempt <= 2:
             try:
                 if attempt == 2:
                     # on second attempt: use pymol to download cif
                     fetch_and_save(id, input[0])
                 parser = MMCIFParser()
                 structure = parser.get_structure(id, input[0])
                 parsing_error = False
             except Exception as e:
                 print(e)
                 attempt += 1
             else:
                 for chain in list(structure.get_models())[0]:
                     header = f">{id.upper()}_{chain.get_id()}"
                     seq = ""
                     for res in chain:
                         res_name = seq1(res.get_resname())
                         seq += res_name if res_name != "X" else ""
                     with open(output[0], "a") as f:
                         fasta_handle = f"{header}\n{seq}"
                         SeqIO.write(SeqIO.parse(StringIO(fasta_handle), "fasta"), f, "fasta")

         if attempt > 2:
             shell("touch {output[0]}")

rule collect:
    input:
         expand("data/temp/fastas_from_cifs/{id}.fasta",
                id=[re.findall(".*?/(\w+)\.cif", path)[0]
                    for path in glob("peptidereactor/db/cifs/*.cif")])
    output:
         "peptidereactor/db/pdb/in_structure/pdb.fasta"
    run:
         with open(output[0], "a") as f_out:
            for path in list(input):
                with open(path) as f_in:
                    records = SeqIO.parse(f_in, "fasta")
                    SeqIO.write(records, f_out, "fasta")

rule make_db:
    input:
         "peptidereactor/db/pdb/in_structure/pdb.fasta"
    output:
         "peptidereactor/db/pdb/in_structure/pdb.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids -blastdb_version 5
         touch {output}
         """