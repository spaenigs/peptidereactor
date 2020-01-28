import sys
sys.path.append("/home/spaenigs/Apps/anaconda3/envs/encoding_benchmark/lib/python3.7/site-packages")

import pandas as pd
import numpy as np
from more_itertools import windowed, tail
from Bio.PDB import PDBParser
from gridData import Grid
from pymol.cgo import *
from pymol import cmd

cmd.load("/home/spaenigs/PycharmProjects/eb/data/temp/096e1693/HAntifreeze.pqr")
cmd.load("/home/spaenigs/PycharmProjects/eb/data/temp/096e1693/HAntifreeze.esp.dx")
points = pd.read_csv("/home/spaenigs/PycharmProjects/eb/data/temp/096e1693/HAntifreeze_6.eh.csv")

# structure = PDBParser()\
#     .get_structure("test", "FULL_PATH_TO_PDB")

g = Grid("/home/spaenigs/PycharmProjects/eb/data/temp/096e1693/HAntifreeze.esp.dx")
sphere_list = []
for x in range(0, 225, 5):
    for y in range(0, 225, 5):
        for z in range(0, 225, 5):
            i, j, k = g.edges[0][x], g.edges[1][y], g.edges[2][z]
            if np.abs(i)+j+k not in range(200):
                 continue
            color = [COLOR, 1, 0, 0] if g.grid[x,y,z] > 0 else [COLOR, 0, 0, 1]
            sphere_list += color
            sphere_list += [SPHERE, i, j, k, 1.5]
            cmd.load_cgo(sphere_list, f'p_{i}_{j}_{k}', 1)

# total_len, required_window_len = len(points["x"]), len(list(structure.get_residues())) - int("ORIGINAL_WINDOW_SIZE")
# res = []
# for ws in [ws for ws in range(1, total_len) if len(str(ws)) < len(str(total_len))]:
# 	for s in [s for s in range(1, total_len) if s < ws]:
# 		windows = enumerate(windowed(range(total_len), ws, step=s), start=1)
# 		size, last_window = list(tail(1, windows))[0]
# 		if last_window[-1] is not None and size == required_window_len:
# 			# 	res += [(total_len, ws, s)]
# 			print((total_len, ws, s))
#
# print(res)

sphere_list = []
for _, (i, j, k) in points.iterrows():
    sphere_list += [COLOR, 255 / 255, 228 / 255, 181 / 255]
    sphere_list += [SPHERE, i, j, k, 0.5]
cmd.load_cgo(sphere_list, f'segment_full', 1)
