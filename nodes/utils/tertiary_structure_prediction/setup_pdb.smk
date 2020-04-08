"""
url_all_ids = \
            "https://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*" + \
            "&customReportColumns=structureId,sequence&format=csv&service=wsfile"

response = requests.get(url)
ids_handle = response.content.decode()
df_ids = pd.read_csv(StringIO(ids_handle), dtype={"structureId": str})
df_ids["structureId"] = df_ids["structureId"].apply(lambda id: id.lower())
df_ids.loc[pd.isna(df_ids["sequence"]) != True, :]


"""

# TODO add rsync and gunzip env yaml

rule download_db:
    output:
         "peptidereactor/db/pdb/pdb.fasta.gz"
    priority:
         50
    threads:
         1000
    run:
         url = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
         shell(f"""
         wget {url} \
            -O {{output}} \
            -q --show-progress --progress=bar:force:noscroll;
         """)

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