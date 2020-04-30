from Bio.Blast.Applications \
    import NcbiblastpCommandline
from modlamp.core import save_fasta, read_fasta
from io import StringIO

import pandas as pd

from peptidereactor.workflow_executer \
    import WorkflowExecuter
from nodes.utils.tertiary_structure_search.scripts.utils \
    import get_seq_names, Cif

include:
    "setup_pdb.smk"

TOKEN = "tyui" # config["token"]
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
         f"data/temp/{TOKEN}/{{seq_name}}.fasta"
    run:
         seqs, names = read_fasta(str(input[0]))
         seq_tuples = dict((name, seq) for name, seq in zip(names, seqs))
         seq = seq_tuples[wildcards.seq_name]
         save_fasta(str(output), sequences=[seq], names=[wildcards.seq_name])

rule blast_search:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.fasta",
         "peptidereactor/db/pdb/in_structure/pdb.db"
    output:
         f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv"
    priority:
        1000
    run:
         header = ["qseqid", "sacc", "sstart", "send", "evalue", "qseq", "sseq"]
         cline = NcbiblastpCommandline(
             task="blastp-short",
             db=str(input[1]),
             query=str(input[0]),
             evalue=200000, # include distinct hits
             max_target_seqs=1,
             max_hsps=1,
             outfmt="10 " + " ".join(header))
         stdout, stderr = cline()

         df_res = pd.read_csv(StringIO(stdout), names=header)
         blast_hits = df_res.shape[0]

         if blast_hits > 0:

             df_res["sacc"] = \
                 df_res["sacc"].apply(lambda full_id: full_id.split("_not_by_pdb_")[0])

             ids, chains = \
                 zip(*df_res["sacc"].apply(
                     lambda full_id: full_id.split("_")).values.tolist())
             df_res["sacc_id"] = [id.lower() for id in ids]
             df_res["sacc_chain"] = chains

             # remove hits with two-letter chain names (see https://stackoverflow.com/questions/50579608)
             df_res = \
                 df_res.loc[ [False if len(id) > 1 else True for id in df_res["sacc_chain"]]
                           , :]

         df_res.to_csv(str(output))

rule utils_download_cif_files:
    input:
         f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv"
    output:
         f"data/temp/{TOKEN}/cifs_downloaded_for_{{seq_name}}.txt"
    params:
         snakefile="nodes/utils/download_cifs/Snakefile",
         configfile="nodes/utils/download_cifs/config.yaml"
    threads:
        1000
    run:
         df = pd.read_csv(str(input[0]), dtype={"sacc_id": str})
         if not df.empty:
             cif_files = expand(CIFS_DIR + "{id}.cif", id=df["sacc_id"])
             with WorkflowExecuter(dict(), dict(cifs_out=cif_files), params.configfile):
                 shell(f"""snakemake -s {{params.snakefile}} --configfile {{params.configfile}} --quiet""")
         shell("touch {output}")

rule dump_cleaved_structure:
    input:
         f"data/temp/{TOKEN}/blast_result_{{seq_name}}.csv",
         f"data/temp/{TOKEN}/cifs_downloaded_for_{{seq_name}}.txt"
    output:
         TARGET_DIR + "{seq_name}.pdb"
    run:
         df = pd.read_csv(str(input[0]), dtype={"sacc_id":str})

         if df.empty:
             shell("touch {output}")
         else:
             pdb_id, chain_id, hit = \
                 df.loc[0, ["sacc_id", "sacc_chain", "sseq"]]

             cif = Cif(pdb_id, chain_id, CIFS_DIR)

             if cif.check_invalid():
                 cif.remove_invalid(f"data/temp/{TOKEN}/{pdb_id}_{chain_id}.pdb")
                 if not cif.invalids_removed:
                     shell("touch {output}")

             cif.dump_slice(hit, output[0])


