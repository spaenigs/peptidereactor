from modlamp.core import read_fasta, save_fasta
from more_itertools import windowed
import joblib as jl
import urllib.parse
import urllib.request
import yaml

TOKEN = config["token"]

def get_ids():
    with open(config["ids_file_in"]) as f:
        return  set([l.rstrip() for l in f.readlines()])

wildcard_constraints:
    id="\w{6}"

rule all:
    input:
         config["fasta_out"],
         config["classes_out"],
         config["fasta_complete_out"],
         config["classes_yaml_out"],
         config["classes_idx_out"]

rule dump_fasta_orig:
    input:
         config["dataset_in"]
    output:
         temp(f"data/temp/{TOKEN}/dataset.jl"),
         config["fasta_out"],
         config["classes_out"]
    run:
         with open(str(input)) as f:
             tuples = [l.rstrip().split(",") for l in f.readlines()]
             seq_tups = dict((kmer, int(class_)) for kmer, class_ in tuples)

         jl.dump(seq_tups, str(output[0]))

         def write_and_flush(stream, s):
             stream.write(s)
             stream.flush()

         with open(str(output[1]), mode="a") as f1, \
                 open(str(output[2]), mode="a") as f2:
             for i, (seq, class_) in enumerate(seq_tups.items(), start=1):
                 write_and_flush(f1, f">Seq_{str(i)}\n{seq}\n")
                 write_and_flush(f2, f"{0 if class_ == -1 else 1}\n")

rule get_fasta:
    input:
         config["ids_file_in"]
    output:
         temp(f"data/temp/{TOKEN}/proteins.fasta")
    run:
         with open(str(input)) as f:
            ids = [l.rstrip() for l in f.readlines()]

         url = "https://www.uniprot.org/uploadlists/"
         params = {"from": "ACC+ID", "to": "ACC", "format": "fasta", "query": " ".join(ids)}
         data = urllib.parse.urlencode(params)
         data = data.encode("utf-8")
         req = urllib.request.Request(url, data)

         with urllib.request.urlopen(req) as f1, \
                 open(str(output), mode="w") as f2:
             response = f1.read().decode("utf-8")
             f2.write(response)
             f2.flush()

rule replace_seq_names:
    input:
         config["ids_file_in"],
         f"data/temp/{TOKEN}/proteins.fasta"
    output:
         config["fasta_complete_out"]
    run:
         with open(str(input[0])) as f:
            ids = [l.rstrip() for l in f.readlines()]

         seqs, names = read_fasta(str(input[1]))

         for i, name in enumerate(names):
             for id in ids:
                 names[i] = id if id in name else names[i]

         save_fasta(str(output), seqs, names)

rule get_windows:
    input:
         config["fasta_complete_out"],
         f"data/temp/{TOKEN}/dataset.jl"
    output:
         config["classes_yaml_out"]
    run:
         seqs, names = read_fasta(str(input[0]))
         res = []
         tmp = []
         for seq, name in zip(seqs, names):

            windows = ["".join(w) for w in windowed(seq, 8)]
            seq_tups = jl.load(str(input[1]))
            kmers = set(set(windows)).intersection(seq_tups.keys())

            tmp_res = {name: []}
            for start, end in zip(range(0, len(seq)+1), range(8, len(seq)+1)):
                window = seq[start:end]
                if window in kmers:
                    # avoid duplicated entries
                    if window in tmp:
                        continue
                    tmp += [window]
                    kmers.remove(window)
                    tmp_res[name] += [{"class": 0 if seq_tups[window] == -1 else 1,
                                       "range": [start, end]}]
            res += [tmp_res]

         with open(str(output), mode="w") as f:
             yaml.safe_dump(res, f)

rule combine_classes_yaml:
    input:
         config["fasta_complete_out"],
         config["classes_yaml_out"]
    output:
         config["classes_idx_out"]
    run:
         seqs, names = read_fasta(str(input[0]))

         with open(input[1]) as f:
            res =yaml.safe_load(f)

         with open(str(output), mode="w") as f:
             d = dict([(list(e.keys())[0], i) for i, e in enumerate(res)])
             for n in names:
                 f.write(f"{d[n]}\n")
                 f.flush()
