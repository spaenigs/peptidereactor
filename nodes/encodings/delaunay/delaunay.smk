from modlamp.core import read_fasta
import pandas as pd
import os

TOKEN = config["token"]
TARGET_FILES = config["csv_out"]
PDB_DIR = config["pdb_dir"]

if type(TARGET_FILES) == list:
    TARGET_DIR = os.path.dirname(TARGET_FILES[0])
else:
    TARGET_DIR = os.path.dirname(TARGET_FILES)

rule all:
    input:
         config["csv_out"]

rule encode:
    input:
         config["fasta_in"]
    output:
         temp(f"data/temp/{TOKEN}/delaunay_{{algorithm}}.csv")
    run:
         from nodes.encodings.delaunay.scripts.delaunay_triangulation import ProteinStructure

         _, names = read_fasta(input[0])

         res = {}
         names_for_index = []
         for name in names:
             pdb_file = PDB_DIR + f"{name}.pdb"
             if wildcards.algorithm == "average_distance":
                 vals = ProteinStructure(protein=pdb_file)\
                     .average_distance()
             elif wildcards.algorithm == "total_distance":
                 vals = ProteinStructure(protein=pdb_file)\
                     .total_distance()
             elif wildcards.algorithm == "number_instances":
                 vals = ProteinStructure(protein=pdb_file)\
                     .number_instances()
             elif wildcards.algorithm == "frequency_instances":
                 vals = ProteinStructure(protein=pdb_file)\
                     .frequency_instances()
             elif wildcards.algorithm == "cartesian_product":
                 vals = ProteinStructure(protein=pdb_file)\
                     .cartesian_product()
             else:
                 raise ValueError(f"Unknown algorithm: {wildcards.algorithm}.")
             if len(vals) > 0:
                 res[name] = vals
                 names_for_index += [name]

         df = pd.DataFrame(res).transpose()
         df.index = names_for_index
         df.to_csv(output[0])

rule dump:
    input:
         f"data/temp/{TOKEN}/delaunay_{{algorithm}}.csv",
         config["fasta_in"],
         config["classes_in"]
    output:
         f"{TARGET_DIR}/delaunay_{{algorithm}}.csv"
    run:
         df = pd.read_csv(input[0], index_col=0)

         seqs, names = read_fasta(input[1])
         with open(input[2]) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         for (name, (seq, class_)) in seq_tuples.items():
             df.loc[name, "y"] = class_

         df.dropna(inplace=True)
         df.to_csv(output[0])
