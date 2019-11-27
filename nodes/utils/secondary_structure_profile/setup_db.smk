import os

# TARGET_DIR = config["target_dir"]

# rule all:
#     input:
#         "apps/db/uniref90/uniref90.db"
#         # TARGET_DIR + "uniref90.db"

rule download_db:
    input:
         config["uniprot90_download_link_in"]
    output:
         "apps/db/uniref90/uniref90.fasta.gz"
         # TARGET_DIR + "uniref90.fasta.gz"
    priority:
         50
    threads:
         1000
    shell:
         """
         wget $(head -n 1 {input[0]}) \
            -P {output} \
            -q --show-progress --progress=bar:force:noscroll;
         """

rule unzip_db:
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
