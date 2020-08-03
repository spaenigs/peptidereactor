from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.stats import gaussian_kde
from itertools import *
from more_itertools import flatten
import numpy as np
import pandas as pd

FUNCTIONAL_TYPES = {
    "N": ["hydrogen_bond_acceptor", "hydrogen_bond_donor"],
    "D": ["hydrogen_bond_acceptor"],
    "Q": ["hydrogen_bond_acceptor", "hydrogen_bond_donor"],
    "E": ["hydrogen_bond_acceptor"],
    "H": ["ambivalent_donor_acceptor", "pi_stacking_centers"],
    "S": ["ambivalent_donor_acceptor"],
    "T": ["ambivalent_donor_acceptor", "aliphatic"],
    "Y": ["ambivalent_donor_acceptor", "pi_stacking_centers"],
    "A": ["aliphatic"],
    "R": ["aliphatic", "hydrogen_bond_donor"],
    "C": ["aliphatic"],
    "I": ["aliphatic"],
    "L": ["aliphatic"],
    "K": ["aliphatic", "hydrogen_bond_donor"],
    "M": ["aliphatic"],
    "P": ["aliphatic"],
    "V": ["aliphatic"],
    "F": ["pi_stacking_centers"],
    "W": ["pi_stacking_centers", "hydrogen_bond_donor"]
}

functional_types = set(flatten(FUNCTIONAL_TYPES.values()))

func_types_tuples = \
    dict((k, v) for k, v in
         zip_longest(functional_types, [], fillvalue=[]))

distances = \
    dict((k, v) for k, v in \
         zip_longest(functional_types, [], fillvalue=func_types_tuples))

parser = PDBParser()
structure = parser.get_structure(snakemake.wildcards.seq_name, snakemake.input[0])

for r1, r2 in combinations(structure.get_residues(), 2):

    aa1, aa2 = seq1(r1.get_resname()), seq1(r2.get_resname())

    if "G" in [aa1, aa2]:
        continue

    if "U" in [aa1, aa2]:
        continue

    for func_type_1, func_type_2 in product(FUNCTIONAL_TYPES[aa1], FUNCTIONAL_TYPES[aa2]):

        def get_carbon_coord(residue):
            # c-alpha
            try:
                return [atom.get_coord() for atom in filter(lambda a: a.get_id() == "CA", residue)][0]
            except IndexError:
                # use index of first atom, if CA not present
                return [atom.get_coord() for atom in residue][0]

        def get_carbon_coord_al(residue):
            coords = [atom.get_coord() for atom in
                      filter(lambda a: "C" in a.get_id() and len(a.get_id()) >= 2, residue)]
            # if no C in residue, take all other atoms
            if len(coords) == 0:
                coords = [atom.get_coord() for atom in residue]
            # aromatic: unweighted average of the respective carbons
            # see http://www.biology.arizona.edu/biochemistry/problem_sets/aa/Aliphatic.html
            return np.array([np.sum(tr) / len(coords) for tr in zip(*coords)])

        def dist_fn(r1, r2):
            return np.linalg.norm(fn1(r1) - fn2(r2))

        if (func_type_1, func_type_2) == ("aliphatic", "aliphatic"):
            fn1, fn2 = get_carbon_coord_al, get_carbon_coord_al
        else:
            if func_type_1 == "aliphatic":
                fn1, fn2 = get_carbon_coord_al, get_carbon_coord
            elif func_type_2 == "aliphatic":
                fn1, fn2 = get_carbon_coord, get_carbon_coord_al
            else:
                fn1, fn2 = get_carbon_coord, get_carbon_coord

        distances[func_type_1][func_type_2] += [dist_fn(r1, r2)]

## compute_distance_distribution:

range_vals = np.arange(0.5, 50.5, 0.5)
res_vec = []
for func_type_1, func_type_2 in combinations_with_replacement(functional_types, 2):
    dataset = distances[func_type_1][func_type_2]
    if len(dataset) <= 1:
        # results to smaller, final dataset
        pass
    elif len(list(set(dataset))) == 1:
        pass
    else:
        kde = gaussian_kde(dataset, bw_method=1)
        res_vec += list(kde.evaluate(range_vals))

pd.DataFrame({snakemake.wildcards.seq_name: res_vec}) \
    .transpose().to_csv(str(snakemake.output))
