rule paac_filter_datasets:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,paac}/csv/filtered/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}.csv"
    group:
        "filter_and_normalize"
    script:
        "../scripts/filter.py"


rule paac_normalize:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/filtered/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,paac}/csv/normalized/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}_normalized-{normalized}.csv"
    group:
        "filter_and_normalize"
    script:
          "../scripts/normalize.py"


def collect_files(wildcards):
    return expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
                        "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}_normalized-{normalized}.csv",
                  dataset=wildcards.dataset, part=wildcards.part,  normalized=wildcards.normalized,
                  encoding=wildcards.encoding,
                  lambdaValue=config["lambda_based"]["paac"]["lambdas"])

rule paac_collect_normalized:
    input:
         collect_files
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,paac}/csv/normalized/" + \
        "{dataset}_{part}_normalized-{normalized}.txt"
    run:
        for path in list(input):
            with open(str(output), mode="a") as f:
                f.write(f"{path}\n")
