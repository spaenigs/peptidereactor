import unittest
import secrets
from snakemake.shell import shell
import pandas as pd
from modlamp.core import read_fasta

DATASET = "test"
TOKEN = secrets.token_hex(4)
BASE_DIR = f"tests/data/{DATASET}"


class TestClassAssignments(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        seqs, names = read_fasta(f"tests/data/{DATASET}/seqs.fasta")
        with open(f"tests/data/{DATASET}/classes.txt") as f:
            classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
        cls.seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))

    def test_class_assignments(self):
        target_file = f"{BASE_DIR}/csv/zscale.csv"
        fasta_in, classes_in = f"{BASE_DIR}/seqs.fasta", f"{BASE_DIR}/classes.txt"
        shell(f"""
        ./apps/run_pipeline -s nodes/encodings/zscale/Snakefile \
            --config dataset={DATASET} token={TOKEN} {target_file} \
                     fasta_in={fasta_in} classes_in={classes_in} csv_out={target_file} \
            --quiet
        """)
        df = pd.read_csv(f"{BASE_DIR}/csv/zscale.csv", index_col=0)
        actual, expected = [], []
        for (name, (seq, class_)) in self.seq_tuples.items():
            actual += [df.loc[name, "y"]]
            expected += [class_]
        self.assertSequenceEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()