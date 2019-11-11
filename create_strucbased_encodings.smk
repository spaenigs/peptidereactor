import os

config["global_workdir"] = os.getcwd() + "/"

DATASET = config["dataset"]

rule encoding_delaunay:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         profile=f"data/{DATASET}/profile"
    output:
         csv_out=expand(f"data/{DATASET}/csv/delaunay/delaunay_{{algorithm}}.csv",
                        algorithm=["average_distance", "total_distance", "cartesian_product",
                                   "number_instances", "frequency_instances"])
    params:
         subworkflow="delaunay",
         snakefile="nodes/encodings/delaunay/Snakefile",
         configfile="nodes/encodings/delaunay/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"