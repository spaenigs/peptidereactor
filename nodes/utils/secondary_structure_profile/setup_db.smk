rule download_db:
    input:
         config["uniprot90_download_link_in"]
    output:
         "apps/db/uniref90/uniref90.fasta.gz"
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
         "apps/db/uniref90/uniref90.fasta.gz"
    output:
         "apps/db/uniref90/uniref90.fasta"
    shell:
         "gunzip -c {input} > {output}"

rule make_db:
    input:
         ancient("apps/db/uniref90/uniref90.fasta")
    output:
         "apps/db/uniref90/uniref90.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids
         touch {output}
         """
