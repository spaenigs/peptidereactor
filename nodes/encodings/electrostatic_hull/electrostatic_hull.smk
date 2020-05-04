from modlamp.core import read_fasta
import pandas as pd
import numpy as np
import yaml
import os

TOKEN = config["token"]
PDB_DIR = config["pdb_dir"]
TARGET_FILES = config["csv_out"]

if type(TARGET_FILES) == list:
    TARGET_DIR = os.path.dirname(TARGET_FILES[0])
else:
    TARGET_DIR = os.path.dirname(TARGET_FILES)

rule all:
    input:
         config["csv_out"]

rule assign_charges_and_radii:
    input:
         PDB_DIR + "{seq_name}.pdb"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.pqr")
    run:
         # if wildcards.seq_name == "Seq_277":
         #     import pydevd_pycharm
         #     pydevd_pycharm.settrace('localhost', port=8889, stdoutToServer=True, stderrToServer=True)
         # shell("pdb2pqr --whitespace --ff=amber --assign-only {input} {output} 1> /dev/null;")
         cmd = "pdb2pqr --whitespace XXX --ff=amber {input[0]} {output[0]} 1> /dev/null;"
         try:
             shell(cmd.replace("XXX", ""))
         except Exception as e:
             try:
                 # Only assign charges and radii - do not add atoms, debump, or optimize
                 shell(cmd.replace("XXX", "--assign-only"))
             except Exception as e2:
                 # Do no optimization, atom addition, or parameter assignment, just return the original
                 # PDB file in aligned format. Overrides --ff
                 shell(cmd.replace("XXX", "--clean"))

rule solvent_accessible_surface:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.pqr"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.sas.dx")
    threads:
         1000
    shell:
         f"nodes/encodings/electrostatic_hull/scripts/run_apbs.sh {TOKEN} {{input}} {{output}} smol"

rule electrostatic_hull:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.sas.dx"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}.eh.csv")
    run:
         from nodes.encodings.electrostatic_hull.scripts.parse_grid \
             import readDX, electrostatic_hull
         dx = readDX(str(input))
         electrostatic_hull(data=dx.data,
                            iterations=int(int(wildcards.distance)/dx.delta[0])+1,
                            filename=str(output))

rule electrostatic_potential:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.pqr"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.esp.dx")
    threads:
         1000
    shell:
         f"nodes/encodings/electrostatic_hull/scripts/run_apbs.sh {TOKEN} {{input}} {{output}} pot"

rule electrostatic_pot_at_electrostatic_hull_grid:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}.eh.csv",
         f"data/temp/{TOKEN}/{{seq_name}}.esp.dx"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_part.csv")
    run:
         from nodes.encodings.electrostatic_hull.scripts.parse_grid \
            import csv2points, readDX, dx2csv

         filter = csv2points(str(input[0]))
         ids, dx_list = [], []
         dx_list.append(readDX(str(input[1])))
         ids.append(wildcards.seq_name)

         dx2csv(dx_list, filter=filter, ids=ids, filename=str(output), sep=",")

rule combine:
    input:
         config["fasta_in"],
         lambda wildcards: \
             expand(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_part.csv",
                    seq_name=read_fasta(config["fasta_in"])[1],
                    distance=wildcards.distance)
    output:
         temp(f"data/temp/{TOKEN}/final_{{distance}}.yaml")
    run:
         enco = {"enco_seqs": {}}
         for path in list(input[1:]):
             df_tmp = pd.read_csv(path, index_col=0)
             seq_name = df_tmp.index[0]
             enco["enco_seqs"][seq_name] = df_tmp.loc[seq_name, :].values.tolist()

         with open(str(output), mode="w") as f:
            enco["interpolate_to"] = \
                int(np.median([len(seq) for seq in read_fasta(str(input[0]))[0]]))
            yaml.safe_dump(enco, f)

rule interpolate:
    input:
         enco=f"data/temp/{TOKEN}/final_{{distance}}.yaml"
    output:
         temp(f"data/temp/{TOKEN}/out_{{distance}}.csv")
    script:
         "scripts/interpolate.R"

rule dump:
    input:
         f"data/temp/{TOKEN}/out_{{distance}}.csv",
         config["fasta_in"],
         config["classes_in"]
    output:
         f"{TARGET_DIR}/electrostatic_hull_{{distance}}.csv"
    run:
         df = pd.read_csv(str(input[0]), index_col=0)
         df["y"] = -1

         seqs, names = read_fasta(str(input[1]))
         fastas = [[n, s] for s, n in zip(seqs, names)]
         with open(str(input[2])) as f:
            classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         for (name, (seq, class_)) in seq_tuples.items():
             df.loc[name, "y"] = class_

         df.sort_values(by="y").to_csv(str(output))
