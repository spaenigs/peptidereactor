import jinja2 as j2

TOKEN = config["token"]

DATASET = config["dataset"]

rule compute_dendrogram:
    input:
         config["metrics_dir_in"] + "f1.csv",
         config["dataset_correlation_in"]
    output:
         config["html_dir_out"] + f"ds_corr/dataset_correlation.json"
    script:
         "scripts/compute_dendrogram.R"

rule create_dendrogram:
    input:
         config["html_dir_out"] + f"ds_corr/dataset_correlation.json"
    output:
         temp(f"data/temp/{TOKEN}/datasets_dendrogram.json")
    run:
         env = j2.Environment(
             loader=j2.FileSystemLoader("nodes/vis/single_dataset/templates/"),
             autoescape=j2.select_autoescape(["json"])
         )

         url = \
             DATASET + "/" + input[0].replace(config["html_dir_out"], "")

         title = [
             "Correlation of encoded datasets, based on the adjusted RV coefficient."
         ]

         template = env.get_template("dendrogram.json")
         template\
            .stream(url=url, title=title)\
            .dump(output[0])
