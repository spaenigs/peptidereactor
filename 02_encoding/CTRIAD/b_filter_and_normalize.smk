rule ctriad_filter_datasets:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_gap-{gapValue}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,ctriad}/csv/filtered/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_gap-{gapValue}.csv"
    group:
        "filter_and_normalize"
    script:
        "../scripts/filter.py"


rule ctriad_normalize:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/filtered/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_gap-{gapValue}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,ctriad}/csv/normalized/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_gap-{gapValue}_normalized-{normalized}.csv"
    group:
        "filter_and_normalize"
    script:
          "../scripts/normalize.py"


def collect_files(wildcards):
    return expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
                        "{dataset}_{part}_ifeature_{encoding}encoder_gap-{gapValue}_normalized-{normalized}.csv",
                  dataset=wildcards.dataset, part=wildcards.part,  normalized=wildcards.normalized,
                  encoding=wildcards.encoding,
                  gapValue=config["gap_based"]["ctriad"]["gaps"])

rule ctriad_collect_normalized:
    input:
         collect_files
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,ctriad}/csv/normalized/{dataset}_{part}_normalized-{normalized}.txt"
    run:
        for path in list(input):
            with open(str(output), mode="a") as f:
                f.write(f"{path}\n")
