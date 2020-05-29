from modlamp.core import read_fasta
from functools import reduce
from sklearn.decomposition import PCA
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import joblib as jl
import pandas as pd
import numpy as np

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
         seqs, names = read_fasta(input[0])
         with open(input[1]) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         seq_tuple = seq_tuples[wildcards.seq_name]
         jl.dump(value=([[wildcards.seq_name, seq_tuple[0]]], seq_tuple[1]),
                 filename=output[0])

# rule find_energy_minimized_conformation:
#     input:
#          PDB_DIR + "{seq_name}.pdb"
#     output:
#          temp(f"data/temp/{TOKEN}/{{seq_name}}_minimized.pdb")
#     threads:
#          int(np.maximum(os.cpu_count() / 8, os.cpu_count()))
#     run:
#          pdb = PDBFile(str(input))
#          forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
#          modeller = Modeller(pdb.topology, pdb.positions)
#          modeller.addHydrogens(forcefield)
#          system = forcefield.createSystem(modeller.topology,
#                                           nonbondedMethod=NoCutoff,
#                                           nonbondedCutoff=1*nanometer,
#                                           constraints=HBonds)
#          integrator = LangevinIntegrator(300*kelvin, 1/picometer, 0.002*picometer)
#          simulation = Simulation(modeller.topology, system, integrator)
#          simulation.context.setPositions(modeller.positions)
#          simulation.minimizeEnergy()
#          simulation.step(10000)
#          positions = simulation.context.getState(getPositions=True).getPositions()
#          PDBFile.writeFile(simulation.topology, positions, open(str(output), 'w'))

rule compute_molecular_descriptors:
    input:
         PDB_DIR + "{seq_name}.pdb",
         # f"data/temp/{TOKEN}/{{seq_name}}_minimized.pdb",
         f"data/temp/{TOKEN}/{{seq_name}}.joblib"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_molecular_descriptors.csv")
    script:
         "scripts/compute_qsar.py"

rule combine:
     input:
          expand(f"data/temp/{TOKEN}/{{seq_name}}_molecular_descriptors.csv",
                seq_name=read_fasta(config["fasta_in"])[1])
     output:
          temp(f"data/temp/{TOKEN}/qsar_raw_combined.csv")
     run:
          # remove class from columns
          unique_colnames = \
              reduce(lambda l1, l2: l1 & l2,
                     [pd.read_csv(path, index_col=0).columns[:-1] for path in list(input)])

          unique_colnames = unique_colnames.append(pd.Index(["y"]))

          df_res = pd.DataFrame()
          for csv_path in list(input):
              df_tmp = pd.read_csv(csv_path, index_col=0)
              df_tmp = df_tmp.loc[:, unique_colnames]
              df_res = pd.concat([df_res, df_tmp])

          df_res.to_csv(output[0])

rule remove_zero_variance:
    input:
         f"data/temp/{TOKEN}/qsar_raw_combined.csv"
    output:
         temp(f"data/temp/{TOKEN}/qsar_raw_removed_zero_variance.csv")
    run:
         df = pd.read_csv(input[0], index_col=0)
         df = df.loc[:, [i != 0.0 for i in df.apply(np.var, axis=0)]]
         df.to_csv(output[0])

rule remove_highly_correlated:
    input:
         f"data/temp/{TOKEN}/qsar_raw_removed_zero_variance.csv"
    output:
         temp(f"data/temp/{TOKEN}/qsar_raw_removed_highly_correlated.csv")
    run:
         df = pd.read_csv(input[0], index_col=0)
         cm = df.corr()

         cmindex, cmcolumns = cm.index.values.tolist(), cm.columns.values.tolist()
         for idx in cmindex:
              for col in cmcolumns:
                  if idx != col and np.abs(cm.loc[idx, col]) >= 0.95:
                      cmindex.remove(col)
                      cmcolumns.remove(col)
                      cm = cm.loc[cmindex, cmcolumns]

         df = df.loc[:, cm.columns.append(pd.Index(["y"]))]
         df.to_csv(output[0])

rule run_pca:
    input:
         f"data/temp/{TOKEN}/qsar_raw_removed_highly_correlated.csv"
    output:
         config["csv_out"]
    run:
         df = pd.read_csv(input[0], index_col=0)
         X, y = df.iloc[:, :-1].values, df["y"]

         pca = PCA(n_components=np.min([7, np.min([X.shape[0], X.shape[1]])]))
         X_new = pca.fit_transform(X, y)

         df_res = pd.DataFrame(X_new, index=df.index)
         df_res["y"] = df["y"]
         df_res.to_csv(output[0])