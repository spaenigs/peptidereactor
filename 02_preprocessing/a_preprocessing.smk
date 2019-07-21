from snakemake.io import temp, expand

from typing import List
from Bio import SeqIO
from sklearn.externals import joblib as jl


localrules: read_fasta, read_classes,
            create_input_data,
            create_normal_distributed_input_data,
            save_as_fasta,
            plot_sequence_length_distribution


rule read_fasta:
    input:
         "01_data/in/fasta/{dataset}.fasta"
    output:
        temp("01_data/out/tmp/{dataset}_from_fasta.joblib")
    run:
        raw_data: List[List[str]] = []
        seq_names = []
        for record in SeqIO.parse(str(input), "fasta"):
            seq_names.append(record.name)
            raw_data.append([record.name, str(record.seq)])
        if len(set(seq_names)) != len(raw_data):
            raise Exception("Seqence names should be unique!")
        jl.dump(value=raw_data, filename=str(output))


rule read_classes:
    input:
        "01_data/in/class/{dataset}_classes.txt"
    output:
        temp("01_data/out/tmp/{dataset}_from_classes.joblib")
    run:
        f = open(str(input), mode="r")
        target = [int(line.rstrip()) for line in f.readlines()]
        jl.dump(value=target, filename=str(output))


rule create_input_data:
    input:
        "01_data/out/tmp/{dataset}_from_fasta.joblib",
        "01_data/out/tmp/{dataset}_from_classes.joblib"
    output:
        temp("01_data/out/tmp/{dataset,[A-Za-z]+}.joblib")
    run:
        from_fasta = jl.load(filename=str(input[0]))
        from_classes = jl.load(filename=str(input[1]))
        sorted_in_da = [tup for tup in sorted(zip(from_fasta, from_classes), key=lambda tup: tup[0][0])]
        from_fasta_sorted = list(map(lambda tup: tup[0], sorted_in_da))
        from_classes_sorted = list(map(lambda tup: tup[1], sorted_in_da))
        input_data = (from_fasta_sorted, from_classes_sorted)
        jl.dump(value=input_data, filename=str(output))


rule create_normal_distributed_input_data:
    input:
        "01_data/out/tmp/{dataset}.joblib"
    output:
        temp("01_data/out/tmp/{dataset,[A-Za-z]+}_normal_distributed.joblib"),
        "01_data/out/{dataset,[A-Za-z]+}/{dataset}_ds1/joblib/{dataset}_ds1_normal_distributed.joblib",
        "01_data/out/{dataset,[A-Za-z]+}/{dataset}_ds2/joblib/{dataset}_ds2_normal_distributed.joblib",
    script:
        "scripts/create_normal_distributed_input_data.py"


checkpoint save_as_fasta:
    input:
         "01_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed.joblib"
    output:
        "01_data/out/{dataset,[A-Za-z]+}/fasta/{dataset}_{part}.fasta",
        "01_data/out/{dataset,[A-Za-z]+}/class/{dataset}_{part}_classes.txt"
    run:
        from modlamp.core import save_fasta
        def create_fasta(input_data_):
            seq_tups = input_data_[0]
            classes = input_data_[1]
            seq_names = list(map(lambda tup: tup[0], seq_tups))
            seqs = list(map(lambda tup: tup[1], seq_tups))
            return seq_names, seqs, classes
        input_data = jl.load(str(input))
        n, s, c = create_fasta(input_data)
        save_fasta(str(output[0]), sequences=s, names=n)
        with open(str(output[1]), mode="w") as f:
            for c_ in c:
                f.write(str(c_) + "\n")
                f.flush()


rule plot_sequence_length_distribution:
    input:
        "01_data/out/tmp/{dataset}.joblib",
        "01_data/out/tmp/{dataset}_normal_distributed.joblib",
         lambda wildcards: expand("01_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed.joblib",
                                  dataset=wildcards.dataset, part=["ds1", "ds2"])
    output:
        "01_data/out/{dataset,[A-Za-z]+}/plots/{dataset}_length_distribution.svg"
    script:
        "scripts/plot_sequence_length_distribution.py"
