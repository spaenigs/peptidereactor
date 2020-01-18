from modlamp.core import read_fasta
from more_itertools import windowed
import joblib as jl
import pandas as pd
import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from rdkit import Chem
from mordred import Calculator, descriptors
from mordred import DetourMatrix, get_descriptors_in_module
from Bio.PDB import PDBParser, Dice

# pipeline adapted from https://www.nature.com/articles/s41598-018-19669-4.pdf

TOKEN = config["token"]
PDB_DIR = config["pdb_dir"]

rule all:
    input:
         config["csv_out"]

rule split_input_data:
    input:
         config["fasta_in"],
         config["classes_in"]
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

rule find_energy_minimized_conformation:
    input:
         PDB_DIR + "{seq_name}.pdb"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_minimized.pdb")
    run:
         from sys import stdout

         pdb = PDBFile(str(input))
         forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
         modeller = Modeller(pdb.topology, pdb.positions)
         modeller.addHydrogens(forcefield)
         system = forcefield.createSystem(modeller.topology,
                                          nonbondedMethod=NoCutoff,
                                          nonbondedCutoff=1*nanometer,
                                          constraints=HBonds)
         integrator = LangevinIntegrator(300*kelvin, 1/picometer, 0.002*picometer)
         simulation = Simulation(modeller.topology, system, integrator)
         simulation.context.setPositions(modeller.positions)
         simulation.minimizeEnergy()
         simulation.reporters.append(StateDataReporter(
             stdout, 1000, step=True, potentialEnergy=True, temperature=True))
         simulation.step(10000)
         positions = simulation.context.getState(getPositions=True).getPositions()
         PDBFile.writeFile(simulation.topology, positions, open(str(output), 'w'))

rule compute_molecular_descriptors:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_minimized.pdb",
         f"data/temp/{TOKEN}/{{seq_name}}.joblib"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_molecular_descriptors.csv")
    run:
         structure = PDBParser()\
             .get_structure(wildcards.seq_name, str(input[0]))
         seqs, classes = jl.load(str(input[1]))

         window_size = 20
         len_residues = len(list(structure.get_residues()))
         df_res = pd.DataFrame()
         for start, end in [[0, len_residues]] \
                 if len_residues <= window_size \
                 else enumerate(range(window_size, len_residues)):
             filename = f"data/temp/{TOKEN}/{structure.get_id()}_{start}_{end}.pdb"
             Dice.extract(structure, "A", start, end, filename)
             pdb = Chem.MolFromPDBFile(filename)
             calc = Calculator(descriptors)
             res = pd.DataFrame(calc.pandas([pdb], quiet=True))
             df_res = pd.concat([df_res, res])

         # compute mean descriptor values over all windows
         df_res = pd.DataFrame(df_res.apply(np.mean)).transpose()
         df_res.index = [wildcards.seq_name]
         # TODO check if classes is correctly assigned
         df_res["y"] = [classes]
         df_res.to_csv(str(output))

rule combine:
     input:
         expand(f"data/temp/{TOKEN}/{{seq_name}}_molecular_descriptors.csv",
                seq_name=read_fasta(config["fasta_in"])[1])
     output:
         temp(f"data/temp/{TOKEN}/qsar_raw.csv")
     run:
         df_res = pd.DataFrame()

         for csv_path in list(input):
             df_tmp = pd.read_csv(csv_path, index_col=0)
             df_res = pd.concat([df_res, df_tmp])

         df_res.to_csv(str(output))

rule remove_invalid:
    input:
         f"data/temp/{TOKEN}/qsar_raw.csv"
    output:
         temp(f"data/temp/{TOKEN}/qsar_raw_removed_invalid.csv")
    run:
         df = pd.read_csv(str(input), index_col=0)
         # watch out for missing values in not all rows, i.e., peptides/proteins
         df = df.loc[:, [s.dtype in [np.float, np.int] for n,s in df.items()]]
         df.dropna(axis="columns", inplace=True)
         df.to_csv(str(output))

rule remove_zero_variance:
    input:
         f"data/temp/{TOKEN}/qsar_raw_removed_invalid.csv"
    output:
         temp(f"data/temp/{TOKEN}/qsar_raw_removed_zero_variance.csv")
    run:
         df = pd.read_csv(str(input), index_col=0)
         df = df.loc[:, [i != 0.0 for i in df.apply(np.var, axis=0)]]
         df.to_csv(str(output))

rule remove_highly_correlated:
    input:
         f"data/temp/{TOKEN}/qsar_raw_removed_zero_variance.csv"
    output:
         temp(f"data/temp/{TOKEN}/qsar_raw_removed_highly_correlated.csv")
    run:
         df = pd.read_csv(str(input), index_col=0)
         cm = df.corr()

         cmindex, cmcolumns = cm.index.values.tolist(), cm.columns.values.tolist()
         for idx in cmindex:
              for col in cmcolumns:
                  if idx != col and np.abs(cm.loc[idx, col]) >= 0.95:
                      cmindex.remove(col)
                      cmcolumns.remove(col)
                      cm = cm.loc[cmindex, cmcolumns]

         df.loc[:, cm.columns.append(pd.Index(["y"]))].to_csv(str(output))

rule run_pca:
    input:
         f"data/temp/{TOKEN}/qsar_raw_removed_highly_correlated.csv"
    output:
         config["csv_out"]
    run:
         from sklearn.decomposition import PCA

         df = pd.read_csv(str(input), index_col=0)
         X, y = df.iloc[:, :-1].values, df["y"]

         pca = PCA(n_components=np.min([7,
                                        np.min([X.shape[0], X.shape[1]])]
                                       ))
         X_new = pca.fit_transform(X, y)

         df_res = pd.DataFrame(X_new, index=df.index)
         df_res["y"] = df["y"]
         df_res.to_csv(str(output))












