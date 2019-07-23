rule psekraac_filter_datasets:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,psekraac}/csv/filtered/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}.csv"
    group:
        "filter_and_normalize"
    script:
        "../scripts/filter.py"


rule psekraac_normalize:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/filtered/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}.csv"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,psekraac}/csv/normalized/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}_normalized-{normalized}.csv"
    group:
        "filter_and_normalize"
    script:
          "../scripts/normalize.py"


def collect_files(wildcards):
    files = []
    for type_ in config["psekraac"]["types"]:
        files += expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
                            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}_normalized-{normalized}.csv",
                        dataset=wildcards.dataset, part=wildcards.part,  normalized=wildcards.normalized,
                        encoding=wildcards.encoding,
                        name=config["psekraac"][type_]["name"],
                        subtype=config["psekraac"][type_]["subtypes"],
                        raactype=config["psekraac"][type_]["raactypes"],
                        ktuple=config["psekraac"][type_]["ktuples"],
                        glambda=config["psekraac"][type_]["glambdas"])
    return files

rule psekraac_collect_normalized:
    input:
         collect_files
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,psekraac}/csv/normalized/{dataset}_{part}_normalized-{normalized}.txt"
    run:
        for path in list(input):
            with open(str(output), mode="a") as f:
                f.write(f"{path}\n")
                f.flush()