import warnings
from io import StringIO
from glob import glob
from Bio import SeqIO, BiopythonWarning
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

import re

warnings.simplefilter('ignore', BiopythonWarning)

rule all:
    input:
         "in_structure.fasta"

rule parse_cif_file:
    input:
         "peptidereactor/db/cifs/{id}.cif"
    output:
         "data/temp/fastas_from_cifs/{id}.fasta"
    run:
         parser = MMCIFParser()
         id = wildcards.id

         try:
             structure = parser.get_structure(id, input[0])

             for chain in list(structure.get_models())[0]:
                 header = f">{id.upper()}_{chain.get_id()}"
                 seq = ""
                 for res in chain:
                     res_name = seq1(res.get_resname())
                     seq += res_name if res_name != "X" else ""
                 with open(output[0], "a") as f:
                     fasta_handle = f"{header}\n{seq}"
                     SeqIO.write(SeqIO.parse(StringIO(fasta_handle), "fasta"), f, "fasta")

         except KeyError as ke:
             print(f"Key {ke} missing in {id}")
             shell("touch {output}")

rule collect:
    input:
         expand("data/temp/fastas_from_cifs/{id}.fasta",
                id=[re.findall(".*?/(\w+)\.cif", path)[0]
                    for path in glob("peptidereactor/db/cifs/*.cif")])
    output:
         "in_structure.fasta"
    run:
         with open(output[0], "a") as f_out:
            for path in list(input):
                with open(path) as f_in:
                    records = SeqIO.parse(f_in, "fasta")
                    SeqIO.write(records, f_out, "fasta")

