import pandas as pd

import time

rule all:
    input:
         expand("c_{dataset}.yaml", dataset=range(3))
         # expand("benchmark_{dataset}.txt", dataset=range(3))

rule myrule_b:
    input:
         "a_{dataset}.yaml"
    output:
         "b_{dataset}.yaml"
    benchmark:
         "myrule_b_{dataset}.benchmark.txt"
    run:
         time.sleep(6)
         shell("touch {output}")

rule myrule_a:
    input:
         "peptidereactor.smk"
    output:
         "a_{dataset}.yaml"
    benchmark:
         "myrule_a_{dataset}.benchmark.txt"
    run:
         time.sleep(12)
         shell("touch {output}")

rule myrule_c:
    benchmark:
         "myrule_c_{dataset}.benchmark.txt"
    input:
         "b_{dataset}.yaml"
    output:
         "c_{dataset}.yaml"
    run:
         time.sleep(8)
         shell("touch {output}")

rule collect_benchmark:
    input:
         lambda wildcards: expand(f"{{name}}_{wildcards.dataset}.benchmark.txt", name=["myrule_a", "myrule_b", "myrule_c"])
    output:
         "benchmark_{dataset}.txt"
    run:
         names = ["myrule_a", "myrule_b", "myrule_c"]

         df_res = pd.DataFrame()
         for name in names:
            df = pd.read_csv(f"{name}_{wildcards.dataset}.benchmark.txt", sep="\t", parse_dates=["h:m:s"])
            df.index = [name]
            df_res = pd.concat([df_res, df])

         df_res.to_csv(output[0])
