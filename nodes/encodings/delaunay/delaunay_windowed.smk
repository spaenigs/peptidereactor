from modlamp.core import read_fasta
from Bio.PDB import PDBParser, Dice
import pandas as pd
import os
import yaml
import joblib as jl

TOKEN = config["token"]
TARGET_FILES = config["csv_out"]
PDB_DIR = config["pdb_dir"]

def run_delaunay(algorithm, pdb_file):
    from nodes.encodings.delaunay.scripts.delaunay_triangulation \
        import ProteinStructure
    if algorithm == "average_distance":
        return ProteinStructure(protein=pdb_file).average_distance()
    elif algorithm == "total_distance":
        return ProteinStructure(protein=pdb_file).total_distance()
    elif algorithm == "number_instances":
        return ProteinStructure(protein=pdb_file).number_instances()
    elif algorithm == "frequency_instances":
        return ProteinStructure(protein=pdb_file).frequency_instances()
    elif algorithm == "cartesian_product":
        return ProteinStructure(protein=pdb_file).cartesian_product()
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}.")

if type(TARGET_FILES) == list:
    TARGET_DIR = os.path.dirname(TARGET_FILES[0])
else:
    TARGET_DIR = os.path.dirname(TARGET_FILES)

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
         temp(f"data/temp/{TOKEN}/delaunay_{{algorithm}}_{{seq_name}}.csv")
    run:
         structure = PDBParser()\
             .get_structure(wildcards.seq_name, str(input[0]))
         chain_id = [c.get_id() for c in structure.get_chains()][0]
         seqs, class_idx = jl.load(str(input[1]))

         with open(str(input[2])) as f:
             windowed_classes = yaml.safe_load(f)

         values = list(windowed_classes[class_idx].values())[0]
         res, names, classes = {}, [], []
         for i, v in enumerate(values, start=1):
             start, end = v["range"]
             name, class_ = f"{structure.get_id()}_part_{str(i)}", v["class"]
             filename = f"data/temp/{TOKEN}/{name}.pdb"
             Dice.extract(structure, chain_id, start, end, filename)
             res[name] = run_delaunay(wildcards.algorithm, filename)
             names += [name]
             classes += [class_]

         df_res = pd.DataFrame(res).transpose()
         df_res.index = names
         df_res["y"] = classes
         df_res.to_csv(str(output))

rule dump:
    input:
         lambda wildcards: \
             expand(f"data/temp/{TOKEN}/delaunay_{wildcards.algorithm}_{{seq_name}}.csv",
                    seq_name=read_fasta(config["fasta_in"])[1])
    output:
         f"{TARGET_DIR}/delaunay_{{algorithm}}.csv"
    run:
         df_res = pd.DataFrame()
         for csv_path in list(input):
             df_res = pd.concat([df_res, pd.read_csv(csv_path, index_col=0)])
         df_res.to_csv(str(output))
