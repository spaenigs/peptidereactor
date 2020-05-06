import pandas as pd

import time
from datetime import datetime, timedelta

rule all:
    input:
         expand("c_{id}.yaml", id=range(3)),
         "benchmark.txt"

rule b:
    input:
         "a_{id}.yaml"
    output:
         "b_{id}.yaml"
    benchmark:
         "b_{id}.benchmark.txt"
    run:
         time.sleep(6)
         shell("touch {output}")

rule a:
    input:
         "peptidereactor.yaml"
    output:
         "a_{id}.yaml"
    benchmark:
         "a_{id}.benchmark.txt"
    run:
         time.sleep(12)
         shell("touch {output}")

rule c:
    input:
         "b_{id}.yaml"
    output:
         "c_{id}.yaml"
    benchmark:
         "c_{id}.benchmark.txt"
    run:
         time.sleep(8)
         shell("touch {output}")

rule collect_benchmark:
    input:
         expand("{name}_{id}.benchmark.txt", name=["a", "b", "c"], id=range(3))
    output:
         "benchmark.txt"
    run:
         names, ids = glob_wildcards("{name}_{id}.benchmark.txt")

         df_res = pd.DataFrame()
         for name in set(names):

             df_tmp = pd.DataFrame()
             for id in set(ids):
                df = pd.read_csv(f"{name}_{id}.benchmark.txt", sep="\t", parse_dates=["h:m:s"])
                df_tmp = pd.concat([df_tmp, df])

             df_tmp_2 = \
                 df_tmp.loc[:, ['s', 'max_rss', 'max_vms', 'max_uss', 'max_pss', 'io_in', 'io_out', 'mean_load']]\
                     .apply(sum, 0).to_frame().transpose()

             df_tmp_2["h:m:s"] = datetime.strptime("00:00:00", "%H:%M:%S")
             for d in df_tmp["h:m:s"]:
                 df_tmp_2["h:m:s"] += timedelta(hours=d.hour, minutes=d.minute, seconds=d.second)

             df_tmp_2.index = [name]
             df_res = pd.concat([df_res, df_tmp_2])

         df_res.to_csv(output[0])

         print(df_res)



