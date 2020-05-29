import re
import subprocess

from rdkit import Chem
from mordred import Calculator, descriptors
from Bio.PDB import PDBParser, Dice
import pandas as pd
import numpy as np
import joblib as jl
from snakemake.shell import shell


def get_pdb_chunks(full_pdb, window_size, len_residues, token):
    len_residues += 1

    for chain in full_pdb.get_chains():
        chain_id = chain.get_id()
        start_idx, stop_idx = \
            list(chain.get_residues())[0].get_id()[1], list(chain.get_residues())[-1].get_id()[1]
        if len_residues <= window_size:
            windows = [[start_idx, stop_idx]]
        else:
            tmp = enumerate(range(window_size, len_residues))
            windows = [(i+start_idx, j+start_idx) for i, j in tmp]
        for start, end in windows:
            filename = f"data/temp/{token}/{full_pdb.get_id()}_{start}_{end}.pdb"
            Dice.extract(full_pdb, chain_id, start, end, filename)

            not_valid = True
            tmp_filename = filename

            # in case several atoms exceed explicit valence
            while not_valid:
                cmd = f"from rdkit import Chem; Chem.MolFromPDBFile('{tmp_filename}')"
                process = subprocess.Popen(["python", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                error_msg = stderr.decode()
                if error_msg == "":
                    not_valid = False
                    mol_obj = Chem.MolFromPDBFile(tmp_filename)
                else:
                    hits = re.findall("Explicit valence for atom # (\d+)", error_msg)
                    if len(hits) > 0:
                        line_number = int(hits[0])
                        import secrets
                        id = secrets.token_hex(5)
                        new_filename = filename.replace(".pdb", f".{id}.pdb")
                        with open(tmp_filename) as f1, open(new_filename, "a") as f2:
                            for i, l in enumerate(f1.readlines()):
                                if i == line_number:
                                    continue
                                else:
                                    f2.write(l)
                                    f2.flush()
                        tmp_filename = new_filename
                    else:
                        raise ValueError(f"Error parsing rdkit error message: {error_msg}!")

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
