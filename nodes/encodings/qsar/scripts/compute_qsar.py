from rdkit import Chem
from mordred import Calculator, descriptors
from Bio.PDB import PDBParser, Dice
from snakemake.shell import shell

import pandas as pd
import numpy as np
import joblib as jl

import warnings


def reset_indices(chain):
    for i, residue in enumerate(chain.get_residues(), start=1):
        res_id = list(residue.id)
        res_id[1] += 1e5  # deal with negative indices
        residue.id = tuple(res_id)
    for i, residue in enumerate(chain.get_residues(), start=1):
        res_id = list(residue.id)
        res_id[1] = i
        residue.id = tuple(res_id)
    return chain


def get_pdb_chunks(full_pdb, window_size, len_residues, token):
    len_residues += 1
    for chain in full_pdb.get_chains():
        import pydevd_pycharm
        pydevd_pycharm.settrace('localhost', port=38473, stdoutToServer=True, stderrToServer=True)
        chain = reset_indices(chain)
        chain_id = chain.get_id()
        start_idx, stop_idx = \
            list(chain.get_residues())[0].get_id()[1], \
            list(chain.get_residues())[-1].get_id()[1]
        if len_residues <= window_size:
            windows = [[start_idx, stop_idx]]
        else:
            windows = enumerate(range(window_size, len_residues))
        for start, end in windows:
            filename = f"data/temp/{token}/{full_pdb.get_id()}_{start}_{end}.pdb"
            Dice.extract(full_pdb, chain_id, start, end, filename)
            mol_obj = Chem.MolFromPDBFile(filename)
            shell(f"rm {filename}")
            yield mol_obj


smk_obj = snakemake

name, path, token = \
    smk_obj.wildcards.seq_name, smk_obj.input[0], smk_obj.config["token"]

structure = PDBParser().get_structure(name, path)
seqs, class_ = jl.load(smk_obj.input[1])

len_residues = len(list(structure.get_residues()))
pdbs = get_pdb_chunks(structure, 20, len_residues, token)

try:
    warnings.filterwarnings("ignore")
    calc = Calculator(descriptors)
    df = calc.pandas(list(pdbs), quiet=True)
    df = df.loc[:, [s.dtype in [np.float, np.int] for n, s in df.items()]]
    df.dropna(axis="columns", inplace=True)
    df = pd.DataFrame(df.apply(np.mean)).transpose()
    df.index = [name]
    df["y"] = [class_]
    df.to_csv(str(smk_obj.output))
except Exception as e:
    print(e)
