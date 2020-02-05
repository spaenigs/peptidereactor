from modlamp.core import read_fasta
from more_itertools import windowed, tail
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
          f"data/temp/{TOKEN}/{{seq_name}}.pqr"
    shell:
          "pdb2pqr --whitespace --ff=amber {input} {output} 1> /dev/null;"

rule solvent_accessible_surface:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.pqr"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.sas.dx"
    threads:
         1000
    shell:
         f"nodes/encodings/electrostatic_hull/scripts/run_apbs.sh {TOKEN} {{input}} {{output}} smol"

rule electrostatic_hull:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.sas.dx"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}.eh.csv"
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
         f"data/temp/{TOKEN}/{{seq_name}}.esp.dx"
    threads:
         1000
    shell:
         f"nodes/encodings/electrostatic_hull/scripts/run_apbs.sh {TOKEN} {{input}} {{output}} pot"

rule electrostatic_pot_at_electrostatic_hull_grid:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}.eh.csv",
         f"data/temp/{TOKEN}/{{seq_name}}.esp.dx"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_part.csv"
    run:
         from nodes.encodings.electrostatic_hull.scripts.parse_grid \
            import csv2points, readDX, dx2csv

         filter = csv2points(str(input[0]))
         ids, dx_list = [], []
         dx_list.append(readDX(str(input[1])))
         ids.append(wildcards.seq_name)

         dx2csv(dx_list, filter=filter, ids=ids, filename=str(output), sep=",")

rule get_window_size:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}.eh.csv",
         config["classes_idx_in"],
         config["classes_in"],
         config["fasta_in"],
    output:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_window_size.yaml"
    run:
         points = pd.read_csv(str(input[0]))

         with open(str(input[1])) as f1, \
                open(str(input[2])) as f2:
             classes_idx = \
                list(map(lambda l: int(l.rstrip()), f1.readlines()))
             windowed_classes = \
                yaml.safe_load(f2)

         seqs, names = read_fasta(str(input[3]))

         seq, name, class_idx = \
             [(s, n, cidx) for s, n, cidx in zip(seqs, names, classes_idx)
              if n == wildcards.seq_name][0]

         values = list(windowed_classes[class_idx].values())[0]

         peptide_len = len(seq)
         window_len = values[0]["range"][-1]
         total_len, required_window_len = len(points["x"]), peptide_len - window_len
         res = []
         for ws in [ws for ws in range(1, total_len) if len(str(ws)) < len(str(total_len))]:
         	 for s in [s for s in range(1, total_len) if s < ws]:
                  windows = enumerate(windowed(range(total_len), ws, step=s), start=1)
                  size, last_window = list(tail(1, windows))[0]
                  if last_window[-1] is not None and size == required_window_len:
                       res += [(ws, s)]

         with open(str(output), mode="w") as f:
             tmp = {"name": name, "seq": seq, "values": values}
             index = -2
             if len(res) == 0:
                 raise ValueError(f"No window size found for {wildcards.seq_name}!")
             elif len(res) == 1:
                 index = 0
             elif len(res) == 2:
                 index = 1
             tmp["window_size"] = res[index]
             yaml.safe_dump(tmp, f)

rule get_windows:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_part.csv",
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_window_size.yaml",
    output:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_windowed.csv"
    run:
         # reformat '12.125.10' to (12, 125, 10)
         df = pd.read_csv(str(input[0]), index_col=0)
         df.columns = [tuple([int(i) for i in cn.split(".")]) for cn in df.columns]

         with open(str(input[1])) as f:
             tmp = yaml.safe_load(f)

         seq, name, values = tmp["seq"], tmp["name"], tmp["values"]

         points_x, window_size, step = \
             [x for x, _, _ in df.columns], \
             tmp["window_size"][0] , \
             tmp["window_size"][1]

         df_res, classes_res = pd.DataFrame(), []
         for i, (w, v) in enumerate(zip(windowed(range(len(points_x)), n=window_size, step=step), values), start=1):
             start, end = w[0], w[-1]
             classes_res += [v["class"]]
             # filter based on window range of eh window
             encoded_seq_window = df.iloc[0, w[0]:w[-1]]
             df_tmp = pd.DataFrame({f"{name}_part_{str(i)}": encoded_seq_window.values})\
                 .transpose()
             df_res = pd.concat([df_res, df_tmp])

         df_res["y"] = classes_res
         df_res.to_csv(str(output))

rule collect:
    input:
         lambda wildcards: \
             expand(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_windowed.csv",
                    seq_name=read_fasta(config["fasta_in"])[1],
                    distance=wildcards.distance)
    output:
         f"data/temp/{TOKEN}/final_{{distance}}.yaml"
    run:
         enco = {"enco_seqs": {}}
         column_sizes = []
         for path in list(input):
             df_tmp = pd.read_csv(path, index_col=0).iloc[:, :-1] # remove class column
             column_sizes += [df_tmp.shape[1]]
             for seq_name, series in df_tmp.iterrows():
                enco["enco_seqs"][seq_name] = df_tmp.loc[seq_name, :].values.tolist()

         with open(str(output), mode="w") as f:
            enco["interpolate_to"] = int(np.median(column_sizes))
            yaml.safe_dump(enco, f)

rule interpolate:
    input:
         enco=f"data/temp/{TOKEN}/final_{{distance}}.yaml",
    output:
         f"data/temp/{TOKEN}/out_{{distance}}.csv"
    script:
        "scripts/interpolate.R"

rule dump:
    input:
         f"data/temp/{TOKEN}/out_{{distance}}.csv",
         lambda wildcards: \
             expand(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_windowed.csv",
                    seq_name=read_fasta(config["fasta_in"])[1],
                    distance=wildcards.distance)
    output:
         f"{TARGET_DIR}/electrostatic_hull_{{distance}}.csv"
    run:
         df = pd.read_csv(str(input[0]), index_col=0)
         df["y"] = -1

         for path in list(input[1:]):
             df_tmp = pd.read_csv(path, index_col=0)
             df.loc[df_tmp.index, "y"] = df_tmp.loc[df_tmp.index, "y"]

         df.to_csv(str(output))

rule set_pymol_script:
    input:
         "nodes/encodings/electrostatic_hull/scripts/visualize_esp.py",
         f"data/temp/{TOKEN}/{{seq_name}}.pqr",
         f"data/temp/{TOKEN}/{{seq_name}}.esp.dx",
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}.eh.csv",
         PDB_DIR + "{seq_name}.pdb",
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_window_size.yaml"
    output:
         PDB_DIR.replace("pdb", "misc/pymol") + f"{{seq_name}}_{{distance}}.py"
    run:
         with open(str(input[5])) as f:
             tmp = yaml.safe_load(f)

         # TODO read and set python lib path
         with open(str(input[0])) as f:
             script_text = f.read()\
                 .replace("FULL_PATH_TO_PQR", str(input[1]))\
                 .replace("FULL_PATH_TO_ESP_DX", str(input[2]))\
                 .replace("FULL_PATH_TO_POINTS", str(input[3]))\
                 .replace("FULL_PATH_TO_PDB", str(input[4]))\
                 .replace("ORIGINAL_WINDOW_SIZE", tmp["values"][0]["range"][-1])\
                 .replace("WINDOW_SIZE", tmp["window_size"][0])\
                 .replace("STEP", tmp["window_size"][0])

         with open(str(output), "w") as f:
            f.write(script_text)



