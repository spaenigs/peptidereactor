rule download_db:
    input:
         config["pdb_fasta_download_link_in"]
    output:
         "peptidereactor/db/pdb/pdb.fasta.gz"
    priority:
         50
    threads:
         1000
    shell:
         """
         wget $(head -n 1 {input[0]}) \
            -O {output} \
            -q --show-progress --progress=bar:force:noscroll;
         """

rule unzip_db:
    input:
         "peptidereactor/db/pdb/pdb.fasta.gz"
    output:
         "peptidereactor/db/pdb/pdb.fasta"
    shell:
         "gunzip -c {input} > {output}"

rule make_db:
    input:
         ancient("peptidereactor/db/pdb/pdb.fasta")
    output:
         "peptidereactor/db/pdb/pdb.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids -blastdb_version 5
         touch {output}
         """