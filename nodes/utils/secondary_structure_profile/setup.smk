import os

TARGET_DIR = config["target_dir"]

rule all:
    input:
        "apps/db/uniref90/uniref90.db"
        # TARGET_DIR + "uniref90.db"

rule download:
    output:
         "apps/db/uniref90/uniref90.fasta.gz"
         # TARGET_DIR + "uniref90.fasta.gz"
    shell:
         "wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz -P {output}"

rule unzip:
    input:
         "apps/db/uniref90/uniref90.fasta.gz"
         # TARGET_DIR + "uniref90.fasta.gz"
    output:
         "apps/db/uniref90/uniref90.fasta"
         # TARGET_DIR + "uniref90.fasta"
    shell:
         "gunzip -c {input} > {output}"

rule make_db:
    input:
         ancient("apps/db/uniref90/uniref90.fasta")
         # ancient(TARGET_DIR + "uniref90.fasta")
    output:
         "apps/db/uniref90/uniref90.db"
         # TARGET_DIR + "uniref90.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids
         touch {output}
         """
