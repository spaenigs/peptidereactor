import json

from nodes.vis.single_dataset.scripts.vega_specs.dendrogram \
    import vega_dendrogram

TOKEN = config["token"]

rule all:
    input:
         config["json_out"]

rule compute_dendrogram:
    input:
         config["metric_csv_in"],
         config["dataset_correlation_in"]
    output:
         f"data/temp/{TOKEN}/dataset_correlation.json"
    script:
         "scripts/compute_dendrogram.R"

rule create_dendrogram:
    input:
         f"data/temp/{TOKEN}/dataset_correlation.json"
    output:
         config["json_out"]
    run:
         with open(input[0]) as f:
             v = vega_dendrogram(json.load(f))

         with open(output[0], "w") as f:
            f.write(json.dumps(v))
            f.flush()


