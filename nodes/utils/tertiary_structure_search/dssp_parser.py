from Bio.PDB import make_dssp_dict

dssp = make_dssp_dict("../../../Seq_4.dssp")

header = """
Dummy file to mimic output from VSL2 Predictor of Intrinsically Disordered Regions

Prediction Scores:
========================================
NO.     RES.    PREDICTION      DISORDER
----------------------------------------
"""

for i, (k, v) in enumerate(dssp[0].items(), start=1):
    print(f"{i} -- {k}:\t{v}")
    header += f"{i}\t{v[0]}\t{100.0 if v[1] == '-' else 0.000}\t{'D' if v[1] == '-' else '.'}\n"

footer = "========================================\n"

with open("../../../Seq_4.dis", "w") as f:
    f.write(header + footer)
    f.flush()





mappings = {
    "H": "H", "G": "H", "I": "H",
    "S": "C", "T": "C", "C": "C", "-": "C",
    "E": "E", "B": "E"
}

header = ["#", "index",  "AA",  "SS", "phi1", "psi1", "P_E", "P_C", "P_H", "phi0", "psi0",
          "ASA", "S_pk", "S_SS", "pk_phi", "pk_psi", "pkc_phi", "pkc_ps"]

table = "\t".join(header) + "\n"
for i, (k, v) in enumerate(dssp[0].items(), start=1):
    ss = mappings[v[1]]
    table += f"\t{i}\t{v[0]}\t{ss}\t{v[3]}\t{v[4]}\t0.0\t0.0\t0.0\t0.0\t0.0\t{v[2]}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n"

with open("../../../Seq_4.spXout", "w") as f:
    f.write(table)
    f.flush()


