from shutil import copyfile

import pandas as pd

import glob
import os


def write_empty_file(path):
    with open(path, mode="w") as f:
        f.write("")
        f.flush()


dataset = snakemake.wildcards.dataset
part = snakemake.wildcards.part
encoding = snakemake.wildcards.encoding
normalized = snakemake.wildcards.normalized

df_gm = pd.read_csv(str(input), index_col=0)

if (not df_gm.empty) or (not df_gm.empty):

    with open(str(snakemake.output), mode="a") as f:

        df_filtered = df_gm.loc[df_gm["min_gm"].apply(lambda x: x == True), :]
        src = glob.glob(
            f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/*{df_filtered.name.values[0]}.csv")[0]
        dst = f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"
        f.write(os.path.basename(src))

        f.flush()

    copyfile(src, dst)

else:
    write_empty_file(str(snakemake.output))