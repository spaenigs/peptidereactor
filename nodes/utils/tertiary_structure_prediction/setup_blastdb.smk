from io import StringIO

import requests
from Bio import SeqIO

rule download_fasta:
    output:
         "peptidereactor/db/swiss_port/proteins.fasta"
    run:
         url = \
             "https://www.uniprot.org/uniprot/?query=length%3A[30+TO+300]+AND+" \
             "reviewed%3Ayes&sort=score&format=fasta"

         fasta_handle = requests.get(url).content.decode()

         with open(output[0], "w") as f:
             sequences = SeqIO.parse(StringIO(fasta_handle), "fasta")
             SeqIO.write(sequences, f, "fasta")

rule make_db:
    input:
         "peptidereactor/db/swiss_port/proteins.fasta"
    output:
         "peptidereactor/db/swiss_port/proteins.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids
         touch {output}
         """
