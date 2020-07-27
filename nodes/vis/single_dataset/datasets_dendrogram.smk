import jinja2 as j2

import json

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
         env = j2.Environment(
             loader=j2.FileSystemLoader("nodes/vis/single_dataset/templates/"),
             autoescape=j2.select_autoescape(["json"])
         )

         with open(input[0]) as f:
             template = env.get_template("dendrogram.json")
             template\
                 .stream(values=json.load(f))\
                 .dump(output[0])
