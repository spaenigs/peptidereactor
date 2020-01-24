from modlamp.core import read_fasta
import pandas as pd
from Bio.PDB import PDBParser, Dice
import yaml
import joblib as jl

TOKEN = config["token"]
PDB_DIR = config["pdb_dir"]

rule all:
    input:
         config["csv_out"]

rule split_input_data:
    input:
         config["fasta_in"],
         config["classes_idx_in"]
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.joblib")
    run:
         seqs, names = read_fasta(str(input[0]))
         with open(str(input[1])) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         seq_tuple = seq_tuples[wildcards.seq_name]
         jl.dump(value=([[wildcards.seq_name, seq_tuple[0]]], seq_tuple[1]),
                 filename=str(output))

rule compute_pdb_chunks:
    input:
         PDB_DIR + "{seq_name}.pdb",
         f"data/temp/{TOKEN}/{{seq_name}}.joblib",
         config["classes_in"]
    output:
         temp(f"data/temp/{TOKEN}/distance_distribution_{{seq_name}}.csv")
    run:
         structure = PDBParser()\
             .get_structure(wildcards.seq_name, str(input[0]))
         chain_id = [c.get_id() for c in structure.get_chains()][0]
         seqs, class_idx = jl.load(str(input[1]))

         with open(str(input[2])) as f:
             windowed_classes = yaml.safe_load(f)

         values = list(windowed_classes[class_idx].values())[0]
         df_res = pd.DataFrame()
         for i, v in enumerate(values, start=1):
             start, end = v["range"]
             name, class_ = f"{structure.get_id()}_part_{str(i)}", v["class"]
             filename = f"data/temp/{TOKEN}/{name}.pdb"
             Dice.extract(structure, chain_id, start, end, filename)
             tmp_csv_path = str(output).replace(".csv", f"_part_{str(i)}.csv")
             shell(f"""python nodes/encodings/distance_distribution/scripts/distance_distribution_windowed.py \
                        {name} {filename} {tmp_csv_path}""")
             df_tmp = pd.read_csv(tmp_csv_path, index_col=0)
             df_tmp["y"] = class_
             df_res = pd.concat([df_res, df_tmp])

         df_res.to_csv(str(output))

rule combine_and_dump:
    input:
         expand(f"data/temp/{TOKEN}/distance_distribution_{{seq_name}}.csv",
                seq_name=read_fasta(config["fasta_in"])[1])
    output:
         config["csv_out"]
    run:
         df_res = pd.DataFrame()
         for csv_path in list(input):
             df_res = pd.concat([df_res, pd.read_csv(csv_path, index_col=0)])
         df_res.to_csv(str(output))
