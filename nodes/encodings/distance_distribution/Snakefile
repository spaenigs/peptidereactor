from modlamp.core import read_fasta
import pandas as pd

TOKEN = config["token"]
PDB_DIR = config["pdb_dir"]

rule all:
    input:
         config["csv_out"]

rule compute_distance_distribution:
    input:
         PDB_DIR + "{seq_name}.pdb"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_distance_distribution.csv")
    script:
         "scripts/distance_distribution.py"

rule dump:
    input:
         config["fasta_in"],
         config["classes_in"],
         expand(f"data/temp/{TOKEN}/{{seq_name}}_distance_distribution.csv",
                seq_name=read_fasta(config["fasta_in"])[1])
    output:
         config["csv_out"]
    run:
         seqs, names = read_fasta(input[0])
         fastas = [[n, s] for s, n in zip(seqs, names)]
         with open(input[1]) as f:
            classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         df = pd.DataFrame()
         for path in list(input[2:]):
             df = pd.concat([df, pd.read_csv(path, index_col=0)])
         df["y"] = -1

         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         for (name, (seq, class_)) in seq_tuples.items():
             df.loc[name, "y"] = class_

         # for short sequences with similar amino acid content the algorithm return NA's
         df.dropna(inplace=True)

         df.sort_values(by="y").to_csv(output[0])
