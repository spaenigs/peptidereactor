import os
import re
import pandas as pd
import requests
from io import StringIO

from Bio import SeqIO
from modlamp.core import save_fasta, read_fasta

TOKEN = config["token"]

wildcard_constraints:
    full_id="\w{4}_\w{1,2}",
    type="in|ex"

def get_response(url):
    response = requests.get(url)
    if response.status_code != 200:
        raise IOError(f"Some issues with PDB now. Try again later...\n(URL: {url}")
    return response.content.decode()

def get_current_ids():
    url = "https://www.rcsb.org/pdb/rest/getCurrent"
    ids_handle = get_response(url)
    matches = \
        [re.search("<PDB structureId=\"(\w+)\" />", line)
         for line in ids_handle.split("\n") if "structureId" in line]
    df_res = pd.DataFrame({"structureId": [m.group(1).lower() for m in matches]})
    return df_res["structureId"].values

def get_full_ids():

    path = "data/temp/full_ids.txt"

    if os.path.exists(path):
        with open(path) as f:
            return [l.rstrip() for l in f.readlines()]

    else:
        url = \
            "https://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*" + \
            "&customReportColumns=structureId,chainId,uniprotAcc,sequence&format=csv&service=wsfile"

        ids_handle = get_response(url)

        df_res = pd.read_csv(StringIO(ids_handle))

        # remove enties \wo sequences
        df_res = df_res.loc[pd.isna(df_res["sequence"]) != True, :]

        df_res["structureId"] = df_res["structureId"].apply(lambda id: id.lower())

        # remove dna sequences
        df_res = df_res.loc[ (pd.isna(df_res["uniprotAcc"]) != True)
                           & [False if len(re.findall("[ATGCX]+", s)) == 1 else True for s in df_res["sequence"]]
                           , :]

        current_ids = get_current_ids()
        df_res = df_res.loc[df_res["structureId"].isin(current_ids), :]

        res = [f"{s['structureId']}_{s['chainId']}" for _, s in df_res.iterrows()]

        with open("data/temp/full_ids.txt", "a") as f:
            for r in res:
                f.write(f"{r}\n")
                f.flush()
        return res

rule all:
    input:
        expand("peptidereactor/db/pdb/{type}_structure/pdb.db",
               type=["in", "ex"])

rule download_sequence:
    output:
         f"data/temp/{TOKEN}/{{full_id}}.fasta"
    run:
         id = wildcards.full_id.upper()

         # https://academic.oup.com/nar/article/36/suppl_2/W255/2506064#42648547
         url = \
             f"http://www.bioinformatics.nl/tools/seqatoms/cgi-bin/getseqs?db=pdb_seqatms&id={id}"

         response = requests.get(url)
         fasta_handle = response.content.decode()

         seqs, names = [], []
         for r in SeqIO.parse(StringIO(fasta_handle), "fasta"):
             seqs += [str(r.seq)]
             names += [r.id]

         save_fasta(str(output), seqs, names)

rule group_sequences:
    input:
         f"data/temp/{TOKEN}/{{full_id}}.fasta"
    output:
         f"data/temp/{TOKEN}/{{full_id}}_in_structure.fasta",
         f"data/temp/{TOKEN}/{{full_id}}_ex_structure.fasta"
    run:
         def concat_seq_parts(seq_lst):
             if len(seq_lst) == 0:
                 return [""]
             else:
                 return ["".join(seq_lst)]

         seqs, names = read_fasta(str(input))

         seq_in_structure_lst = re.findall("([A-Z]+)", seqs[0])
         save_fasta(str(output[0]), concat_seq_parts(seq_in_structure_lst), names)

         seq_ex_structure_lst = re.findall("([a-z]+)", seqs[0])
         save_fasta(str(output[1]), concat_seq_parts(seq_ex_structure_lst), names)

rule collect_fastas:
    input:
         lambda wildcards:
             expand(f"data/temp/{TOKEN}/{{full_id}}_{wildcards.type}_structure.fasta",
                    full_id=get_full_ids())
    output:
         "peptidereactor/db/pdb/{type}_structure/pdb.fasta"
    run:
         seqs, names = [], []
         for path in list(input):
             seqs_tmp, names_tmp = read_fasta(path)
             if len(seqs_tmp[0]) == 0:
                 continue
             if names_tmp[0] in names:
                 continue
             seqs += seqs_tmp
             names += names_tmp

         save_fasta(str(output), seqs, names)

rule make_db:
    input:
         "peptidereactor/db/pdb/{type}_structure/pdb.fasta"
    output:
         "peptidereactor/db/pdb/{type}_structure/pdb.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids -blastdb_version 5
         touch {output}
         """