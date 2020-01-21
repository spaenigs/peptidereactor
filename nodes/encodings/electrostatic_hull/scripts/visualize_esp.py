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

import pandas as pd
import numpy as np
from more_itertools import windowed
from pymol.cgo import COLOR, SPHERE
from pymol import cmd

points = pd.read_csv("data/temp/d2f6eb6b/UniRef100_A6P3B2_0.eh.csv")

sphere_list = []
for _, (i, j, k) in points.iterrows():
	sphere_list += [COLOR, 255/255, 228/255, 181/255]
	sphere_list += [SPHERE, i, j, k, 0.5]
cmd.load_cgo(sphere_list, f'segment_full', 1)

w = list(windowed(set(points["x"]), 20))[0]
print(w)
r, g, b = np.random.randint(255)/255, np.random.randint(255)/255, np.random.randint(255)/255
sphere_list = []
for _, (x, y, z) in points.loc[points["x"].isin(range(w[0], w[-1])), :].iterrows():
	sphere_list += [COLOR, 50/255, 205/255, 50/255]
	sphere_list += [SPHERE, x, y, z, 1.00]
cmd.load_cgo(sphere_list, f'segment_{w[0]}_{w[-1]}',   1)
