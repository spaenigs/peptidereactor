from tempfile import NamedTemporaryFile

import pandas as pd
import altair as alt
import numpy as np
import jinja2 as j2

import altair_saver
import re
import json
import joblib

from nodes.vis.single_dataset.scripts.utils import cluster, is_struc_based

TOKEN = config["token"]


def compute_median_f1(path):
    dataset = re.findall("data/(.*?)/", path)[0]
    return pd\
        .read_csv(path, index_col=0)\
        .apply(np.median)\
        .to_frame(dataset)


def is_imbalanced(path):
    with open(path) as f:
        classes = [int(i.rstrip()) for i in f.readlines()]
        return sum(classes) / len(classes)


rule overview_data:
    input:
        config["metrics_dirs_in"]
    output:
        temp(f"data/temp/{TOKEN}/hm_data.csv"),
        temp(f"data/temp/{TOKEN}/hm_dendro_data.csv")
    run:
        paths = [p + "f1.csv" for p in list(input)]

        df_res = pd.DataFrame()
        for p in paths:
            df = pd.read_csv(p, index_col=0)
            df_medians = df.apply(np.median).to_frame("median")
            group = lambda enc: \
                "psekraac" if "lambda-corr" in enc or "g-gap" in enc else enc[:6]
            df_medians["group"] = [group(x) for x in df_medians.index]

            df_tmp = df_medians.groupby(by="group").max()
            df_tmp.columns = [re.findall("data/(.*?)/", p)[0]]
            df_res = pd.concat([df_res, df_tmp], axis=1)

        df_res.fillna(0.0, inplace=True)

        row_indices = cluster(df_res.values, axis=0)
        col_indices = cluster(df_res.values, axis=1)
        heatmap_data = df_res.iloc[row_indices, col_indices]

        df_res.to_csv(output[0])
        heatmap_data.to_csv(output[1])

rule source_data:
    input:
         f"data/temp/{TOKEN}/hm_data.csv"
    output:
         config["html_dir_out"] + f"data/overview_hm/hm_source.json"
    run:
         df_res = pd.read_csv(input[0], index_col=0)

         x, y = np.meshgrid(df_res.columns, df_res.index)

         source = pd.DataFrame({"Dataset": x.ravel(), "Encoding": y.ravel(), "F1": df_res.values.ravel()})

         source["type"] = ["structure based" if is_struc_based(e) else "sequence based" for e in source["Encoding"]]
         source["is_imbalanced"] = source["Dataset"].apply(lambda ds: is_imbalanced(f"data/{ds}/classes.txt"))
         source.sort_values(by="is_imbalanced", inplace=True)

         fields = sorted(source["Dataset"].apply(lambda x: x[:3]).unique())
         fields_dict = dict((k, v) for k, v in zip(fields, range(1, len(fields)+1)))
         source["bio_field"] = source["Dataset"].apply(lambda ds: ds[:3])

         type_dict = {"sequence based": 1, "structure based": 2}
         source["type_field"] = source["type"].apply(lambda t: type_dict[t])

         source.to_json(output[0], orient="records")

rule heatmap_biomedical_application:
    input:
         config["html_dir_out"] + f"data/overview_hm/hm_source.json",
    output:
         temp(f"data/temp/{TOKEN}/hm_bio_appl.joblib")
    run:
         source = pd.read_json(input[0])

         sorted_encodings = source\
             .groupby(by="Encoding")\
             .mean().sort_values("F1", ascending=False)\
             .index.to_list()

         hcharts = []
         for b in sorted(source["bio_field"].unique()):
             vcharts = []
             for t in source["type"].unique():
                 df_tmp = source.loc[source["bio_field"] == b, :].loc[source["type"] == t, :]
                 df_tmp.sort_values(by="is_imbalanced", inplace=True)
                 full_url = input[0].replace("source", f"source_{b}_{t.replace(' ', '')}")
                 df_tmp.to_json(full_url, orient="records")
                 url = input[0].replace(config["html_dir_out"], "").replace("source", f"source_{b}_{t.replace(' ', '')}")
                 vcharts += [alt.Chart(url).mark_rect(strokeWidth=5.0).encode(
                         y=alt.Y(
                             'Encoding:N',
                             title=None if not b.startswith("ace") else \
                                 "Sequence-based encodings" if t == "sequence based" else "Structure-based encodings",
                             axis=alt.Axis(labels=False if not b.startswith("ace") else True,
                                           ticks=False if not b.startswith("ace") else True),
                             sort=alt.Sort(alt.SortArray(sorted(sorted_encodings)))
                         ),
                         x=alt.Y(
                             'Dataset:N',
                             title=None,
                             axis=alt.Axis(labels=False if t == "sequence based" else True,
                                           ticks=False if t == "sequence based" else True),
                             sort=alt.Sort(alt.SortArray(source["Dataset"].unique()))
                         ),
                         color=alt.Color("F1:Q", title="F1"),
                         tooltip=["Encoding:N", "Dataset:N", "F1:Q", "is_imbalanced:Q"]
                     ).properties(
                         height=500 if t == "sequence based" else 150
                     )
                 ]
             hcharts += [alt.vconcat(*vcharts, spacing=2)]

         joblib.dump(alt.hconcat(*hcharts, spacing=5), output[0])

rule heatmap_imbalanced:
    input:
         config["html_dir_out"] + f"data/overview_hm/hm_source.json",
    output:
         temp(f"data/temp/{TOKEN}/hm_imbalanced.joblib")
    run:
         source = pd.read_json(input[0])

         fields_dict = {"sequence based": 1, "structure based": 2}
         source["type_field"] = source["type"].apply(lambda t: fields_dict[t])

         names_imbalanced = source["Dataset"].unique()

         vcharts = []
         for t in source["type"].unique():
             hcharts = []
             for thresholds in [[0.0, 0.35], [0.35, 1.0]]:
                 df_tmp = source.loc[source["type"] == t, :].loc[source["is_imbalanced"].between(*thresholds), :]
                 full_url = input[0].replace("source", f"source_{t}_{thresholds[0]}")
                 df_tmp.to_json(full_url, orient="records")
                 url = input[0].replace(config["html_dir_out"], "").replace("source", f"source_{t}_{thresholds[0]}")
                 hcharts += [alt.Chart(url).mark_rect().encode(
                     y=alt.Y(
                         'Encoding:N',
                         title=None,
                         # title=None if thresholds[1] == 1.0 else \
                         #     "Sequence-based encodings" if t == "sequence based" else "Structure-based encodings",
                         axis=alt.Axis(labels=False, ticks=False)
                         # axis=alt.Axis(labels=False if thresholds[1] == 1.0 else True,
                         #               ticks=False if thresholds[1] == 1.0 else True),
                     ),
                     x=alt.X(
                         'Dataset:N',
                         sort=alt.Sort(alt.SortArray(names_imbalanced)),
                         title=None,
                         # title=None if t == "sequence based" else ["Imbalanced", "datasets"] if thresholds[1] == 0.3 else [
                         #     "Balanced", "datasets"],
                         axis=alt.Axis(labels=False if t == "sequence based" else True,
                                       ticks=False if t == "sequence based" else True)
                     ),
                     color=alt.Color("F1:Q", title="F1"),
                     tooltip=["Encoding:N", "Dataset:N", "F1:Q", "is_imbalanced:Q"]
                 ).properties(
                     height=500 if t == "sequence based" else 150
                 )]
             vcharts += [alt.hconcat(*hcharts, spacing=2)]

         joblib.dump(alt.vconcat(*vcharts, spacing=2), output[0])

rule concat_heatmaps:
    input:
         f"data/temp/{TOKEN}/hm_bio_appl.joblib",
         f"data/temp/{TOKEN}/hm_imbalanced.joblib",
         f"data/temp/{TOKEN}/hm_cluster.json"
    output:
         temp(f"data/temp/{TOKEN}/overview_hm.json")
    run:
         hm_bio = joblib.load(input[0])
         hm_imb = joblib.load(input[1])

         chart = alt.vconcat(hm_bio, hm_imb)

         with open(output[0], "w") as f:
             f.write(chart.to_json(indent=None))
             f.flush()

rule make_dendro_data:
    input:
         f"data/temp/{TOKEN}/hm_dendro_data.csv"
    output:
         config["html_dir_out"] + "data/overview_hm/hm_dendro_{axis}_data.json"
    script:
         "scripts/compute_dendrogram.R"

rule cluster_hm_data:
    input:
         f"data/temp/{TOKEN}/hm_dendro_data.csv"
    output:
         config["html_dir_out"] + f"data/overview_hm/hm_dendro_data_hm_values.json",
         temp(f"data/temp/{TOKEN}/hm_dendro_data_0_transform.json"),
         temp(f"data/temp/{TOKEN}/hm_dendro_data_1_transform.json")
    run:
         heatmap_data = pd.read_csv(input[0], index_col=0)
         x, y = np.meshgrid(heatmap_data.columns, heatmap_data.index)

         source = pd.DataFrame({"Dataset": x.ravel(),
                                "Encoding": y.ravel(),
                                "F1": heatmap_data.values.ravel()
                                })

         source.to_json(output[0], orient="records")
         url = output[0].replace(config["html_dir_out"], "")

         chart = alt.Chart(source).mark_rect().encode(
             x=alt.X("Dataset:N", sort=alt.Sort(alt.SortArray(heatmap_data.columns.to_list()))),
             y=alt.X("Encoding:N", sort=alt.Sort(alt.SortArray(heatmap_data.index.to_list()))),
             color=alt.condition(alt.datum.F1 == 0.0, alt.value("white"), "F1:Q"),
             tooltip="F1:Q"
         )

         with NamedTemporaryFile("w", dir=f"data/temp/{TOKEN}/", suffix=".vg.json") as f:
             altair_saver.save(chart, f.name)
             with open(f.name) as json_in, \
                     open(output[0], "w") as data_hm_out, \
                     open(output[1], "w") as data_0_out, \
                     open(output[2], "w") as data_1_out:
                 vg_spec = json.load(json_in)
                 vg_spec["data"][0]["name"] = "data-hm"
                 vg_spec["data"][1]["source"] = "data-hm"
                 json.dump(vg_spec["data"][0]["values"], data_hm_out),
                 json.dump(vg_spec["data"][1]["transform"], data_0_out),
                 json.dump(vg_spec["data"][2]["transform"], data_1_out)

rule aggregate_hm_cluster:
    input:
         expand(config["html_dir_out"] + "data/overview_hm/hm_dendro_{axis}_data.json",
                axis=["row", "col"]),
         config["html_dir_out"] + f"data/overview_hm/hm_dendro_data_hm_values.json",
         f"data/temp/{TOKEN}/hm_dendro_data_0_transform.json",
         f"data/temp/{TOKEN}/hm_dendro_data_1_transform.json"
    output:
         temp(f"data/temp/{TOKEN}/hm_cluster.json")
    run:
         env = j2.Environment(
             loader=j2.FileSystemLoader("nodes/vis/multiple_datasets/templates/"),
             autoescape=j2.select_autoescape(["json"])
         )

         with open(input[0]) as f:
             data = json.load(f)

         url_row_dendro = input[0].replace(config["html_dir_out"], "")
         url_col_dendro = input[1].replace(config["html_dir_out"], "")
         url_hm_values = input[2].replace(config["html_dir_out"], "")

         with open(input[3]) as f1, open(input[4]) as f2:
             data_0_transform = json.load(f1)
             data_1_transform = json.load(f2)

         template = env.get_template("hm_cluster.json")
         template\
            .stream(url_row_dendro=url_row_dendro,
                    url_col_dendro=url_col_dendro,
                    url_hm_values=url_hm_values,
                    data_0_transform=data_0_transform,
                    data_1_transform=data_1_transform)\
            .dump(output[0])
