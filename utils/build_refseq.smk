import os

TARGET_DIR = config["target_dir"]

rule all:
    input:
        TARGET_DIR + "uniref90.db"

if not os.path.exists(TARGET_DIR + "uniref90.fasta.gz"):
    rule download: 
        output:
            TARGET_DIR + "uniref90.fasta.gz"
        shell:
            f"wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz -P {TARGET_DIR}"

rule unzip:
    input:
         TARGET_DIR + "uniref90.fasta.gz"
    output:
         TARGET_DIR + "uniref90.fasta"
    shell:
         "gunzip -c {input} > {output}"

rule make_db:
    input:
         ancient(TARGET_DIR + "uniref90.fasta")
    output:
         TARGET_DIR + "uniref90.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids
         touch {output}
         """
