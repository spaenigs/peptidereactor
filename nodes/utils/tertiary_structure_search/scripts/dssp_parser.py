from Bio.PDB import make_dssp_dict

import textwrap


def generate_disorder_profile(path_to_dssp, path_to_dis):
    header = textwrap.dedent(
        """\
            Dummy file to mimic output from VSL2 Predictor of Intrinsically Disordered Regions
            
            Prediction Scores:
            ========================================
            NO.     RES.    PREDICTION      DISORDER
            ----------------------------------------
        """)
    with open(path_to_dis, "w") as f:
        dssp = make_dssp_dict(path_to_dssp)
        for i, (k, v) in enumerate(dssp[0].items(), start=1):
            # print(f"{i} -- {k}:\t{v}")
            header += f"{i}\t{v[0]}\t{100.0 if v[1] == '-' else 0.000}\t{'D' if v[1] == '-' else '.'}\n"
        footer = "========================================\n"
        f.write(header + footer)
        f.flush()


def generate_spinex_profile(path_to_dssp, path_to_spXout):
    mappings = {
        "H": "H", "G": "H", "I": "H",
        "S": "C", "T": "C", "C": "C", "-": "C",
        "E": "E", "B": "E"
    }
    header = ["#", "index",  "AA",  "SS", "phi1", "psi1", "P_E", "P_C", "P_H", "phi0", "psi0",
              "ASA", "S_pk", "S_SS", "pk_phi", "pk_psi", "pkc_phi", "pkc_ps"]
    with open(path_to_spXout, "w") as f:
        dssp = make_dssp_dict(path_to_dssp)
        table = "\t".join(header) + "\n"
        for i, (k, v) in enumerate(dssp[0].items(), start=1):
            ss = mappings[v[1]]
            table += f"\t{i}\t{v[0]}\t{ss}\t{v[3]}\t{v[4]}\t0.0\t0.0\t0.0\t0.0\t0.0\t{v[2]}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n"
        f.write(table)
        f.flush()
