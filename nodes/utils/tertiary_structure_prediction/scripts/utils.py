from Bio.PDB import Select


class XDeselect(Select):
    def accept_residue(self, residue):
        if residue.get_resname() == "MSE":
            return 0
        elif residue.get_resname() == "UNK":
            return 0
        else:
            return 1
