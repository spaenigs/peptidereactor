import unittest
import secrets
from snakemake.shell import shell
import pandas as pd
from modlamp.core import read_fasta

TOKEN = secrets.token_hex(4)


class TestClassAssignments(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        seqs, names = read_fasta(f"tests/seqs.fasta")
        with open(f"tests/classes.txt") as f:
            classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
        cls.seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))

    def test_class_assignments(self):
        target_file = "tests/csv/zscale.csv"
        fasta_in, classes_in = "tests/seqs.fasta", "tests/classes.txt"
        shell(f"""
        ./apps/run_pipeline -s nodes/encodings/zscale/Snakefile {target_file} \
            --config token={TOKEN} fasta_in={fasta_in} classes_in={classes_in} csv_out={target_file} \
            --quiet
        """)
        df = pd.read_csv("tests/csv/zscale.csv", index_col=0)
        actual, expected = [], []
        for (name, (seq, class_)) in self.seq_tuples.items():
            actual += [df.loc[name, "y"]]
            expected += [class_]
        self.assertSequenceEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()