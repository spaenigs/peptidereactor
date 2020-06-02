from rdkit import Chem
from mordred import Calculator, descriptors
from Bio.PDB import PDBParser, Dice
from snakemake.shell import shell

import pandas as pd
import numpy as np
import joblib as jl


def get_pdb_chunks(full_pdb, window_size, len_residues, token):
    len_residues += 1
    for chain in full_pdb.get_chains():
        chain_id = chain.get_id()
        start_idx, stop_idx = \
            list(chain.get_residues())[0].get_id()[1], \
            list(chain.get_residues())[-1].get_id()[1]
        if len_residues <= window_size:
            windows = [[start_idx, stop_idx]]
        else:
            tmp = enumerate(range(window_size, len_residues))
            windows = [(i+start_idx, j+start_idx) for i, j in tmp]
        for start, end in windows:
            filename = f"data/temp/{token}/{full_pdb.get_id()}_{start}_{end}.pdb"
            Dice.extract(full_pdb, chain_id, start, end, filename)
            mol_obj = Chem.MolFromPDBFile(filename)
            yield filename, mol_obj


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

try:
    df = calc.pandas(list(pdbs), quiet=True)
    df = df.loc[:, [s.dtype in [np.float, np.int] for n, s in df.items()]]
    df.dropna(axis="columns", inplace=True)
    df = pd.DataFrame(df.apply(np.mean)).transpose()
    df.index = [name]
    df["y"] = [class_]

    df.to_csv(str(smk_obj.output))

    clean_up(filenames)
except Exception as e:
    print(e)
