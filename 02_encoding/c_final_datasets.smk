import sys
sys.path.append("02_encoding")

import scripts.utils as utils

rule generate_distance_matrix:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
        "{dataset}_{part}_{type}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/correlation/" + \
        "{dataset}_{part}_{type}_normalized-{normalized}_vs_rest.csv"
    script:
        "scripts/generate_distance_matrix.py"


rule collect_distance_matrix:
    input:
        lambda wildcards: \
            expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/correlation/" + \
                   "{dataset}_{part}_{type}_normalized-{normalized}_vs_rest.csv",
                   dataset=wildcards.dataset,
                   part=wildcards.part,
                   encoding=wildcards.encoding,
                   normalized=wildcards.normalized,
                   type=utils.get_type(wildcards.encoding, config))
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/" + \
          "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    script:
        "scripts/collect_distance_matrix.py"


rule run_clustering:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/" + \
        "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized,yes|no}.csv"
    run:
        from scripts.run_clustering import run_clustering
        run_clustering(utils.ENCODING_PATTERN[wildcards.encoding],
                       str(input),
                       str(output))


rule compute_geometric_median:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/geom_median/" + \
        "{dataset}_{part}_normalized-{normalized,yes|no}_{type}_vs_rest.csv"
    script:
        "scripts/compute_geometric_median.py"


rule collect_geometric_median:
    input:
        lambda wildcards: \
            expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/geom_median/" + \
                   "{dataset}_{part}_normalized-{normalized}_{type}_vs_rest.csv",
                   dataset=wildcards.dataset,
                   part=wildcards.part,
                   normalized=wildcards.normalized,
                   encoding=wildcards.encoding,
                   type=utils.get_unique_types(wildcards.encoding))
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
         "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    script:
        "scripts/collect_geometric_median.py"


rule plot_clustering:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        "00_data/out/{dataset}/plots/{dataset}_{part}_{encoding}_normalized-{normalized}_tsne.svg"
    script:
        "scripts/plot_clustering.py"


rule get_final_datasets:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        temp("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" +
             "geom_median/tsne/normalized-{normalized}/final_datasets.txt")
    script:
         "scripts/get_final_datasets.py"