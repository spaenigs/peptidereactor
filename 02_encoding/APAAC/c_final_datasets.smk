from snakemake.io import expand, temp

import sys

sys.path.append("02_encoding/")


rule generate_distance_matrix:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}_normalized-{normalized}.csv",
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
        "{dataset}_{part}_normalized-{normalized}.txt"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/correlation/" + \
            "{dataset}_ifeature_{encoding}encoder_lambda-{lambdaValue}_normalized-{normalized}_vs_rest.csv"
    script:
        "../scripts/generate_distance_matrix.py"


def collect_files(wildcards):
    return expand( "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/correlation/" + \
                   "{dataset}_ifeature_{encoding}encoder_lambda-{lambdaValue}_normalized-{normalized}_vs_rest.csv",
                  dataset=wildcards.dataset, part=wildcards.part,  normalized=wildcards.normalized,
                  encoding=wildcards.encoding,
                  lambdaValue=config["lambda_based"]["apaac"]["lambdas"])

rule collect_distance_matrix:
    input:
        collect_files
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/" + \
          "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    script:
         "../scripts/collect_distance_matrix.py"


rule run_clustering:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/" + \
        "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized,yes|no}.csv"
    run:
       from scripts.run_clustering import run_clustering
       run_clustering("apaacencoder_(.*)", str(input), str(output))


rule compute_geometric_median:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/geom_median/" + \
        "{dataset}_{part}_normalized-{normalized,yes|no}_{type}_vs_rest.csv"
    script:
        "..scripts/compute_geometric_median.py"



rule collect_geometric_median:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/geom_median/" + \
                                 "{dataset}_{part}_normalized-{normalized}_{type}_vs_rest.csv",
                                 dataset=wildcards.dataset, part=wildcards.part,
                                 normalized=wildcards.normalized, encoding=wildcards.encoding,
                                 type=["lambda"])
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
         "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    script:
        "../scripts/collect_geometric_median.py"


rule plot_clustering:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        "00_data/out/{dataset}/plots/{dataset}_{part}_{encoding}_normalized-{normalized}_tsne.svg"
    script:
        "..scripts/plot_clustering.py"


rule get_final_datasets:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        temp("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/geom_median/tsne/" + \
             "normalized-{normalized}/final_datasets.txt")
    script:
          "../scripts/get_final_datasets.py"
