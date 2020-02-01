TOKEN = config["token"]

def get_ids():
    with open(config["ids_file_in"]) as f:
        return  set([l.rstrip() for l in f.readlines()])

wildcard_constraints:
    id="\w{6}"

rule all:
    input:
         config["fasta_out"],
         config["classes_out"]

rule dump_fasta_orig:
    input:
         config["dataset_in"]
    output:
         config["fasta_out"],
         config["classes_out"]
    run:
         with open(str(input)) as f:
             tuples = [l.rstrip().split(",") for l in f.readlines()]
             seq_tups = dict((kmer, int(class_)) for kmer, class_ in tuples)

         def write_and_flush(stream, s):
             stream.write(s)
             stream.flush()

         with open(str(output[0]), mode="a") as f1, \
                 open(str(output[1]), mode="a") as f2:
             for i, (seq, class_) in enumerate(seq_tups.items(), start=1):
                 write_and_flush(f1, f">Seq_{str(i)}\n{seq}\n")
                 write_and_flush(f2, f"{0 if class_ == -1 else 1}\n")
