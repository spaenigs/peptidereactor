TARGET_DIR = config["target_dir"]

rule all:
    input:
        TARGET_DIR + "uniref90.fasta"

rule download:
    output:
         TARGET_DIR + "uniref90.fasta.gz"
    shell:
         "wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"

rule unzip:
    input:
         ancient(TARGET_DIR + "uniref90.fasta.gz")
    output:
         TARGET_DIR + "uniref90.fasta"
    shell:
         # https://superuser.com/questions/139419/how-do-i-gunzip-to-a-different-destination-directory/139422
         "gunzip -c {input} > {output}"


rule make_db:
    input:
         TARGET_DIR + "uniref90.fasta"
    output:
         "",
         ""
    shell:
         # https://www.exoscale.com/syslog/bioinformatics-managing-blast-data/
         "makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids"
