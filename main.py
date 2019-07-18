import subprocess, sys

pre_datasets = ["neuropeptides"]
datasets = ["neuropeptides_ds1"]
encodings = ["psekraac"]

## preprocessing
for ds in pre_datasets:
    subprocess.run(
        ["snakemake", "--snakefile", f"02_preprocessing/preprocessing.sf", "--config", f"dataset={ds}", "--cores", "4"])
    for p in ["ds1"]:
        subprocess.run(
            ["snakemake", "--snakefile", f"02_preprocessing/pssm.sf", "--config", f"dataset={ds}", f"part={p}", "--cores", "4"])

## encoding -> filtering -> machine learning
# for ds in datasets:
#     for e in encodings:
#         subprocess.run(["snakemake", "-S", f"encoding/{e}/{e}.sf", "--config", f"dataset={ds}", "--cores", "4"])


