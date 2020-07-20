import json

from nodes.vis.single_dataset.scripts.vega_specs.dendrogram \
    import vega_dendrogram

TOKEN = config["token"]

rule compute_dendrogram:
    input:
         config["metrics_dir_in"] + "f1.csv",
         config["dataset_correlation_in"]
    output:
         temp(f"data/temp/{TOKEN}/dataset_correlation.json")
    script:
         "scripts/compute_dendrogram.R"

rule create_dendrogram:
    input:
         f"data/temp/{TOKEN}/dataset_correlation.json"
    output:
         temp(f"data/temp/{TOKEN}/datasets_dendrogram.json")
    run:
         with open(input[0]) as f:
             v = vega_dendrogram(json.load(f))

         with open(output[0], "w") as f:
            f.write(json.dumps(v))
            f.flush()


