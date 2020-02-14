from modlamp.core import read_fasta
from more_itertools import windowed, split_before
import pandas as pd
import numpy as np
import yaml
import os
from scipy import interpolate

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
    shell:
          "pdb2pqr --whitespace --ff=amber {input} {output} 1> /dev/null;"

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
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_part.csv"
    run:
         from nodes.encodings.electrostatic_hull.scripts.parse_grid \
            import csv2points, readDX, dx2csv

         filter = csv2points(str(input[0]))
         ids, dx_list = [], []
         dx_list.append(readDX(str(input[1])))
         ids.append(wildcards.seq_name)

         dx2csv(dx_list, filter=filter, ids=ids, filename=str(output), sep=",")

rule scale_sliding_windows_to_electrostatic_hull:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_part.csv",
         config["classes_idx_in"],
         config["classes_in"],
         config["fasta_in"]
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_windowed.csv")
    run:
         def get_ws_s(peptide_len, ws, s, eh_points_len):
             factor = eh_points_len / peptide_len
             if factor < 1:
                 zeros = \
                     split_before(
                         np.format_float_positional(factor).replace(".", ""),
                         lambda digit: int(digit) != 0)
                 factor *= eval(f"10e{sum([1 for digit in list(zeros)[0] if digit == '0']) - 1}")
             factor = int(np.round(factor))
             return peptide_len * factor, ws*factor, s*factor

         def interpolate_eh_points(eh_points_len_old, eh_points, eh_points_len_new):
             x, y = np.arange(0, eh_points_len_old), eh_points
             f = interpolate.interp1d(x, y, kind="linear")
             return f(np.linspace(0, eh_points_len_old-1, eh_points_len_new))

         with open(str(input[1])) as f1, \
                open(str(input[2])) as f2:

             classes_idx, windowed_classes = \
                list(map(lambda l: int(l.rstrip()), f1.readlines())), \
                yaml.safe_load(f2)

             seqs, names = read_fasta(str(input[3]))
             seq, name, class_idx = \
                 [(s, n, cidx) for s, n, cidx in zip(seqs, names, classes_idx)
                  if n == wildcards.seq_name][0]

             values = list(windowed_classes[class_idx].values())[0]

         eh_points = pd.read_csv(str(input[0]), index_col=0)
         ws, s = 8, 1
         peptide_len = len(seq)
         eh_points_len = eh_points.shape[1]

         eh_points_len_new, ws_new, s_new = \
             get_ws_s(peptide_len, ws, s, eh_points_len)
         eh_points_interp = \
             interpolate_eh_points(eh_points_len, eh_points.iloc[0, :].values, eh_points_len_new)

         windows_peptide, windows_peptide_indices = \
             windowed(seq, ws, step=s), \
             windowed(range(peptide_len), ws, step=s)
         windows_eh_points, windows_eh_points_indices = \
             windowed(eh_points_interp, ws_new, step=s_new), \
             windowed(range(eh_points_len_new), ws_new, step=s_new)

         df_res, classes_res, cnt = pd.DataFrame(), [], 1
         for wpi, wp, wehi, weh in \
                 zip(windows_peptide_indices,
                     windows_peptide,
                     windows_eh_points_indices,
                     windows_eh_points):
             for ([start, stop], class_) in [(i["range"], i["class"]) for i in values]:
                 if [wpi[0], wpi[-1]] == [start, stop-1]:
                     classes_res += [class_]
                     df_tmp = pd.DataFrame({f"{name}_part_{str(cnt)}": weh})
                     df_res = pd.concat([df_res, df_tmp.transpose()])
                     cnt += 1

         df_res["y"] = classes_res
         df_res.to_csv(str(output))

rule collect:
    input:
         lambda wildcards: \
             expand(f"data/temp/{TOKEN}/{{seq_name}}_{{distance}}_windowed.csv",
                    seq_name=read_fasta(config["fasta_in"])[1],
                    distance=wildcards.distance)
    output:
         temp(f"data/temp/{TOKEN}/final_{{distance}}.yaml")
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
         temp(f"data/temp/{TOKEN}/out_{{distance}}.csv")
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



