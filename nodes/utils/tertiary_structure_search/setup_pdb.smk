from Bio import SeqIO
from modlamp.core import save_fasta, read_fasta
from more_itertools import chunked
from io import StringIO

import pandas as pd

import datetime
import re

from nodes.utils.tertiary_structure_search.scripts.utils \
    import get_response

TOKEN = config["token"]

wildcard_constraints:
    full_id="\w{4}_\w{1,4}",
    type="in|ex"

rule setup_pdb:
    input:
         expand("peptidereactor/db/pdb/{type}_structure_sa/pdb.db",
                type=["in", "ex"])

rule get_ids:
    output:
         f"data/temp/{TOKEN}/full_ids.txt"
    run:
         def get_current_ids():
             url = "https://www.rcsb.org/pdb/rest/getCurrent"
             ids_handle = get_response(url)
             matches = \
                 [re.search("<PDB structureId=\"(\w+)\" />", line)
                  for line in ids_handle.split("\n") if "structureId" in line]
             df_res = pd.DataFrame({"structureId": [m.group(1) for m in matches]})
             return df_res["structureId"].values

         columns = "structureId,chainId,uniprotAcc,releaseDate,sequence"
         url = \
            "https://www.rcsb.org/pdb/rest/customReport.csv" \
            "?pdbids=*" \
            f"&customReportColumns={columns}" \
            "&format=csv&service=wsfile"

         ids_handle = get_response(url)
         df_res = pd.read_csv(StringIO(ids_handle))

         # remove enties \wo sequences
         df_res = df_res.loc[pd.isna(df_res["sequence"]) != True, :]
         # remove dna sequences
         df_res = df_res.loc[ (pd.isna(df_res["uniprotAcc"]) != True)
                            & [False if len(re.findall("[ATGCX]+", s)) == 1 else True for s in df_res["sequence"]]
                            , :]
         # remove sequences with len < 4
         df_res = df_res.loc[[len(s) > 3 for s in df_res["sequence"]], :]
         # dropping ALL duplicate sequences
         df_res.drop_duplicates(subset ="sequence", inplace = True)
         # remove all sequences which not have been added to the seqatoms db yet
         seqatoms_db_version = datetime.date(2020, 3, 23)
         df_res["releaseDate"] = df_res["releaseDate"]\
             .apply(lambda d: datetime.date(*[int(di) for di in d.split("-")]))
         df_res = df_res.loc[df_res["releaseDate"] < seqatoms_db_version]
         # skip obsolete structures
         current_ids = get_current_ids()
         df_res = df_res.loc[df_res["structureId"].isin(current_ids), :]

         # convert to correct chain id: 'NA'
         df_res.loc[pd.isna(df_res["chainId"]), "chainId"] = "NA"

         with open(str(output), "a") as f:
             res = [f"{s['structureId'].upper()}_{s['chainId']}" for _, s in df_res.iterrows()]
             for r in res:
                 f.write(f"{r}\n")
                 f.flush()

rule download_sequences:
    input:
         f"data/temp/{TOKEN}/full_ids.txt"
    output:
         f"data/temp/{TOKEN}/pdb_masked.fasta"
    run:
         with open(str(input)) as f:
             full_ids = [l.rstrip() for l in f.readlines()]

         def get_seqs(id_list):

             # https://academic.oup.com/nar/article/36/suppl_2/W255/2506064#42648547
             url = \
                f"http://www.bioinformatics.nl/tools/seqatoms/cgi-bin/getseqs?db=pdb_seqatms&id=" + ",".join(id_list)

             fasta_handle = get_response(url)

             if "Software error" in fasta_handle:
                 bad_id = re.findall("Could\snot\sextract\ssequence:\s(\w+)", fasta_handle)[0]
                 print(bad_id)
                 id_list.remove(bad_id)
                 get_seqs(id_list)

             seqs_tmp, names_tmp = [], []
             for r in SeqIO.parse(StringIO(fasta_handle), "fasta"):
                 seqs_tmp += [str(r.seq)]
                 names_tmp += [r.id.replace("pdbsa|", "")]

             return seqs_tmp, names_tmp

         seqs, names = [], []
         for chunk in list(chunked(full_ids, 250)):
             seqs_chunk, names_chunk = get_seqs(chunk)
             seqs += seqs_chunk
             names += names_chunk

         save_fasta(str(output), seqs, names)

rule remove_duplicates:
    input:
         f"data/temp/{TOKEN}/pdb_masked.fasta"
    output:
         "peptidereactor/db/pdb/pdb_masked.fasta"
    run:
         seqs, names = read_fasta(str(input))

         seqs_res, names_res = [], []
         all_names, cnt = [], 1
         for seq, name in zip(seqs, names):
             # makeblastdb is case insensitive: 1FNT_C == 1FNT_c
             if name.upper() in all_names:
                 name += "_not_by_pdb_" + str(cnt)
                 cnt += 1
             if seq not in seqs_res:
                 seqs_res += [seq]
                 names_res += [name]
                 all_names += [name.upper()]

         save_fasta(str(output), seqs_res, names_res)

rule group_sequences:
    input:
         f"peptidereactor/db/pdb/pdb_masked.fasta"
    output:
         "peptidereactor/db/pdb/in_structure_sa/pdb.fasta",
         "peptidereactor/db/pdb/ex_structure_sa/pdb.fasta"
    run:
         def concat_seq_parts(seq_lst):
             if len(seq_lst) == 0:
                 return ""
             else:
                 return "".join(seq_lst)

         seqs, names = read_fasta(str(input))

         seqs_in, names_in = [], []
         seqs_ex, names_ex = [], []
         all_names = []
         cnt = 1
         for seq, name in zip(seqs, names):

             seq_in_structure_lst = re.findall("([A-Z]+)", seq)
             seq_in_tmp = concat_seq_parts(seq_in_structure_lst)
             if len(seq_in_tmp) != 0:
                 seqs_in += [seq_in_tmp]
                 names_in += [name]

             seq_ex_structure_lst = re.findall("([a-z]+)", seq)
             seq_ex_tmp = concat_seq_parts(seq_ex_structure_lst)
             if len(seq_ex_tmp) != 0:
                 seqs_ex += [seq_ex_tmp]
                 names_ex += [name]

             all_names += [name.upper()]

         save_fasta(str(output[0]), seqs_in, names_in)
         save_fasta(str(output[1]), seqs_ex, names_ex)

rule make_db:
    input:
         "peptidereactor/db/pdb/{type}_structure_sa/pdb.fasta"
    output:
         "peptidereactor/db/pdb/{type}_structure_sa/pdb.db"
    shell:
         """
         makeblastdb -dbtype prot -in {input} -out {output} -parse_seqids -blastdb_version 5
         touch {output}
         """