"""
To visualize the electrostatic potential by PyMol.

# open e.g. UniRef100_A6P3B2.pqr
# open e.g. UniRef100_A6P3B2.esp.dx
# run Plugin -> APBS Electrostatics -> Run (requires PyMol license)

PyMol> run visualize_esp.py
PyMol> set seq_view, on
PyMol> cmd.delete("segment*")
PyMol> png test, dpi=300
"""

import sys
sys.path.append("/home/spaenigs/Apps/miniconda3/envs/encoding_benchmark/lib/python3.7/site-packages")

import pandas as pd
import numpy as np
from more_itertools import windowed, tail
from Bio.PDB import PDBParser
from pymol.cgo import *
from pymol import cmd

# cmd.load("/media/spaenigs/4B1DB7375F3291A1/peptidereactor/data/temp/096e1693/HAntifreeze.pqr")
cmd.load("FULL_PATH_TO_PQR")
# cmd.load("/media/spaenigs/4B1DB7375F3291A1/peptidereactor/data/temp/096e1693/HAntifreeze.esp.dx")
cmd.load("FULL_PATH_TO_ESP_DX")
# points = pd.read_csv("/media/spaenigs/4B1DB7375F3291A1/peptidereactor/data/temp/096e1693/HAntifreeze_6.eh.csv")
points = pd.read_csv("FULL_PATH_TO_POINTS")

# w = 0.06 # cylinder width
# l = 0.75 # cylinder length
# h = 0.25 # cone hight
# d = w * 1.618 # cone base diameter

# obj = [
# 	CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
#     CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
#     CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
#     CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
#     CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
#     CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0
# ]
#
# cmd.load_cgo(obj, 'axes')
# cmd.load_cgo([COLOR, 1, 0, 0, SPHERE, -3.070870e+02, -2.976540e+02, -3.140975e+02, 5.5, ], "p", 1)
# cmd.load_cgo([COLOR, 1, 0, 0, SPHERE, -7.162213e-03, -7.202935e-03, -7.243608e-03, 1.5, ], "p2", 1)
# cmd.load_cgo([COLOR, 1, 0, 0, SPHERE, 0, 0, 0, 1.5, ], "p3", 1)

#
# cmd.load_cgo([ SPHERE, 0, 0, 0, 0.5, ], f'origin', 1)
# cmd.load_cgo([COLOR, 1, 0, 0, SPHERE, 10, 0, 0, 0.5, ], f'+x', 1)
# # cmd.load_cgo([COLOR, 1, 0, 0, SPHERE, -10, 0, 0, 0.5, ], f'-x', 1)
#
# cmd.load_cgo([COLOR, 0, 1, 0, SPHERE, 0, 10, 0, 0.5, ], f'+y', 1)
# cmd.load_cgo([COLOR, 0, 1, 0, SPHERE, 0, -10, 0, 0.5, ], f'-y', 1)
#
# cmd.load_cgo([COLOR, 0, 0, 1, SPHERE, 0, 0, 10, 0.5, ], f'+z', 1)
# cmd.load_cgo([COLOR, 0, 0, 1, SPHERE, 0, 0, -10, 0.5, ], f'-z', 1)

structure = PDBParser()\
    .get_structure("test", "FULL_PATH_TO_PDB")

total_len, required_window_len = len(points["x"]), len(list(structure.get_residues())) - int("ORIGINAL_WINDOW_SIZE")
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

# for w in windowed(range(total_len), 65, step=14):
for w in windowed(range(total_len), int("WINDOW_SIZE"), step=int("STEP")):
	r, g, b = np.random.randint(255) / 255, np.random.randint(255) / 255, np.random.randint(255) / 255
	sphere_list = []
	for _, (x, y, z) in points.loc[w[0]:w[-1], :].iterrows():
		sphere_list += [COLOR, 50 / 255, 205 / 255, 50 / 255]
		sphere_list += [SPHERE, x, y - 100, z - 80, 1.00]
	cmd.load_cgo(sphere_list, f'segment_{w[0]}_{w[-1]}', 1)

# w = list(windowed(range(len(points["x"])), 492, 9))[-1]
# print(w)
# r, g, b = np.random.randint(255)/255, np.random.randint(255)/255, np.random.randint(255)/255
# sphere_list = []
# for _, (x, y, z) in points.loc[w[0]:w[-1], :].iterrows():
# 	sphere_list += [COLOR, 50/255, 205/255, 50/255]
# 	sphere_list += [SPHERE, x+400, y+40, z-180, 1.00]
# cmd.load_cgo(sphere_list, f'segment_{w[0]}_{w[-1]}',   1)

# for w in windowed(set(points["x"]), 20):
# 	print(len([i for i in points.loc[points["x"].isin(range(w[0], w[-1])), :].iterrows()]))
