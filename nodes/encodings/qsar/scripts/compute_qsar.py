from rdkit import Chem
from mordred import Calculator, descriptors
from Bio.PDB import PDBParser, Dice
import pandas as pd
import numpy as np
import joblib as jl
from snakemake.shell import shell


def get_pdb_chunks(full_pdb, window_size, len_residues, token):
    len_residues += 1
    for start, end in [[0, len_residues]] \
            if len_residues <= window_size \
            else enumerate(range(window_size, len_residues)):
        filename = f"data/temp/{token}/{full_pdb.get_id()}_{start}_{end}.pdb"
        Dice.extract(full_pdb, "A", start, end, filename)
        yield filename, Chem.MolFromPDBFile(filename)


def clean_up(filenames):
    for path in filenames:
        shell(f"rm {path}")


smk_obj = snakemake

name, path, token = \
    smk_obj.wildcards.seq_name, smk_obj.input[0], smk_obj.config["token"]

structure = PDBParser().get_structure(name, path)
seqs, class_ = jl.load(smk_obj.input[1])

calc = Calculator(descriptors)
len_residues = len(list(structure.get_residues()))
filenames, pdbs = zip(*get_pdb_chunks(structure, 20, len_residues, token))

df = calc.pandas(list(pdbs), quiet=True)
df = df.loc[:, [s.dtype in [np.float, np.int] for n, s in df.items()]]
df.dropna(axis="columns", inplace=True)
df = pd.DataFrame(df.apply(np.mean)).transpose()
df.index = [name]
df["y"] = [class_]

df.to_csv(str(smk_obj.output))

clean_up(filenames)
