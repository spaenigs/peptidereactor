from rdkit import Chem
from mordred import Calculator, descriptors
from Bio.PDB import PDBParser, Dice
import joblib as jl
import numpy as np
import yaml
from snakemake.shell import shell


def get_pdb_chunks(full_pdb, windowed_classes, token):
    values = list(windowed_classes[class_idx].values())[0]
    for i, v in enumerate(values, start=1):
        start, end = v["range"]
        name, class_ = f"{full_pdb.get_id()}_part_{str(i)}", v["class"]
        filename = f"data/temp/{token}/{full_pdb.get_id()}_{start}_{end}.pdb"
        Dice.extract(full_pdb, "A", start, end, filename)
        yield name, class_, filename, Chem.MolFromPDBFile(filename)


def clean_up(filenames):
    for path in filenames:
        shell(f"rm {path}")


smj_obj = snakemake

structure = PDBParser()\
    .get_structure(smj_obj.wildcards.seq_name, str(smj_obj.input[0]))
seqs, class_idx = jl.load(str(smj_obj.input[1]))

with open(str(smj_obj.input[2])) as f:
    windowed_classes = yaml.safe_load(f)

names, classes, filenames, pdbs = \
    zip(*get_pdb_chunks(structure, windowed_classes, smj_obj.config["token"]))

calc = Calculator(descriptors)
df = calc.pandas(list(pdbs), quiet=True)
df = df.loc[:, [s.dtype in [np.float, np.int] for n, s in df.items()]]
df.dropna(axis="columns", inplace=True)
df.index, df["y"] = names, classes
df.to_csv(str(smj_obj.output))

clean_up(filenames)
