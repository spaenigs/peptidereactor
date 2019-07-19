import subprocess, sys

pre_datasets = ["neuropeptides"]
datasets = ["neuropeptides_ds1"]
encodings = ["psekraac"]

## preprocessing
for ds in pre_datasets:

    cmd = ["snakemake", "--snakefile", f"02_preprocessing/preprocessing.sf", "--config", f"dataset={ds}",
           "--cores", "4", "-np"]
    if sys.argv[1] == "":
        cmd.pop()
    subprocess.run(cmd)

    for p in ["ds1"]:
        cmd = ["snakemake", "--snakefile", f"02_preprocessing/pssm.sf", "--config", f"dataset={ds}", f"part={p}",
               "--cores", "8", "-np"]
        if sys.argv[1] == "":
            cmd.pop()
        subprocess.run(cmd)

## encoding -> filtering -> machine learning
# for ds in datasets:
#     for e in encodings:
#         subprocess.run(["snakemake", "-S", f"encoding/{e}/{e}.sf", "--config", f"dataset={ds}", "--cores", "4"])


