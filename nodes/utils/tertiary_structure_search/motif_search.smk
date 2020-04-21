from io import StringIO
from Bio import SeqIO
from modlamp.core \
    import save_fasta, read_fasta

import pandas as pd

import sys
import os

sys.path.append(".")

from peptidereactor.workflow_executer import \
    WorkflowExecuter
from nodes.utils.tertiary_structure_search.scripts.utils \
    import dump_structure_slice, get_seq_names

include:
    "setup_pdb.smk"

TOKEN = config["token"]
CIFS_DIR = "peptidereactor/db/cifs/"
TARGET_DIR = config["pdbs_out"]

rule all:
    input:
        expand(TARGET_DIR + "{seq_name}.pdb",
               seq_name=get_seq_names(config["fasta_in"]))

rule split_input_data:
    input:
         config["fasta_in"]
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.fasta")
    run:
         seqs, names = read_fasta(str(input[0]))
         seq_tuples = dict((name, seq) for name, seq in zip(names, seqs))
         seq = seq_tuples[wildcards.seq_name]
         save_fasta(str(output), sequences=[seq], names=[wildcards.seq_name])

rule motif_search:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         "peptidereactor/db/pdb/in_structure/pdb.fasta",
         "peptidereactor/db/pdb/in_structure/pdb.db"
    output:
         temp(f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv")
    run:
         seqs, names = read_fasta(str(input[0]))
         seq, name = seqs[0], names[0]

         cmd = f"cat {str(input[1])} | grep {seq} -B 1"

         handle = os.popen(cmd)\
             .read()\
             .rstrip()\
             .replace("--", "")

         df_res = pd.DataFrame()
         for r in SeqIO.parse(StringIO(handle), "fasta"):
             pdb_id, chain_id = r.id.split("_not_by_pdb_")[0].split("_")
             if len(chain_id) > 1:
                 continue
             start = str(r.seq).find(seq)
             end = start + 8
             df_tmp = pd.DataFrame({
                 "qseqid": [name], "sacc": [r.id],
                 "sstart": [start], "send": [end], "evalue": [-1],
                 "qseq": [seq], "sseq": [seq],
                 "sacc_id": [pdb_id.lower()], "sacc_chain": [chain_id]})
             df_res = pd.concat([df_res, df_tmp])

         df_res.to_csv(str(output))

rule utils_download_cif_files:
    input:
         f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv"
    output:
         temp(f"data/temp/{TOKEN}/cifs_downloaded_for_{{seq_name}}.txt")
    params:
         snakefile="nodes/utils/download_cifs/Snakefile",
         configfile="nodes/utils/download_cifs/config.yaml"
    threads:
         1000
    run:
         df = pd.read_csv(str(input[0]), dtype={"sacc_id": str})
         if not df.empty:
             cif_files = expand(CIFS_DIR + "{id}.cif", id=df["sacc_id"])
             with WorkflowExecuter(dict(), dict(cifs_out=cif_files), params.configfile) as e:
                 shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
         shell("touch {output}")

rule dump_cleaved_structure:
    input:
         f"data/temp/{TOKEN}/motse_result_{{seq_name}}.csv",
         f"data/temp/{TOKEN}/cifs_downloaded_for_{{seq_name}}.txt"
    output:
         TARGET_DIR + "{seq_name}.pdb"
    run:
         df = pd.read_csv(str(input[0]), dtype={"sacc_id":str, "sacc_chain": str})
         if df.empty:
             shell("touch {output}")
         else:
             pdb_id, chain_id, hit = \
                 df.loc[0, ["sacc_id", "sacc_chain", "sseq"]]
             dump_structure_slice(pdb_id, chain_id, hit, CIFS_DIR, str(output))
