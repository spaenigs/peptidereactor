import scripts.utils as utils

localrules: filter_datasets, normalize_datasets, collect_normalized_datasets

rule filter_datasets:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
          "{dataset}_{part}_{type}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/filtered/" + \
         "{dataset}_{part}_{type}.csv"
    group:
        "filter_and_normalize"
    script:
        "scripts/filter.py"


rule normalize_datasets:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/filtered/" + \
            "{dataset}_{part}_{type}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
            "{dataset}_{part}_{type}.csv"
    group:
        "filter_and_normalize"
    script:
          "scripts/normalize.py"


rule collect_normalized_datasets:
    input:
         lambda wildcards: \
            expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
                  "{dataset}_{part}_{type}.csv",
                  dataset=wildcards.dataset,
                  part=wildcards.part,
                  encoding=wildcards.encoding,
                  type=utils.get_type(wildcards.encoding, config))
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/{dataset}_{part}_normalized-{normalized}.txt"
    run:
        for path in list(input):
            with open(str(output), mode="a") as f:
                f.write(f"{path}\n")
