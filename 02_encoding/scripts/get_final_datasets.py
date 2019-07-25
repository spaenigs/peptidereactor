from shutil import copyfile

import pandas as pd

import glob
import os
import sys

sys.path.append(os.getcwd())
sys.path.append(os.getcwd() + "/apps")

import scripts.utils as utils
from get_aaindex_files import get_aaindex_files

def write_empty_file(path):
    with open(path, mode="w") as f:
        f.write("")
        f.flush()


dataset = snakemake.wildcards.dataset
part = snakemake.wildcards.part
encoding = snakemake.wildcards.encoding
normalized = snakemake.wildcards.normalized

if encoding == utils.AAINDEX:

    df = pd.read_csv(str(snakemake.input), index_col = 0)
    get_aaindex_files(df, str(snakemake.output), dataset, part, normalized)

elif encoding in [utils.APAAC, utils.PAAC, utils.CKSAAGP, utils.CKSAAP, utils.CTRIAD, utils.KSCTRIAD,
                            utils.GEARY, utils.MORAN, utils.NMBROTO, utils.QSORDER, utils.SOCNUMBER,
                            utils.EAAC, utils.EGAAC, utils.PSEKRAAC]:

    df_gm = pd.read_csv(str(snakemake.input), index_col=0)

    if not df_gm.empty:

        with open(str(snakemake.output), mode="a") as f:

            df_filtered = df_gm.loc[df_gm["min_gm"].apply(lambda x: x == True), :]
            src = glob.glob(
                f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
                f"*{df_filtered.name.values[0]}.csv")[0]
            dst = f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
                  f"geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"

            f.write(f"{src},{dst}\n")
            f.flush()

        copyfile(src, dst)

    else:
        write_empty_file(str(snakemake.output))

elif encoding in [utils.BINARY, utils.AAC, utils.GAAC, utils.CTDT, utils.CTDC, utils.CTDD,
                  utils.TPC, utils.DPC, utils.GDPC, utils.DDE, utils.BLOSUM62, utils.ZSCALE,
                  utils.GTPC]:

    with open(str(snakemake.output), mode="a") as f:

        src = glob.glob(
            f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/*.csv")[0]
        dst = f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
              f"geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"

        copyfile(src, dst)

        f.write(f"{src},{dst}\n")
        f.flush()

else:
    raise ValueError(f"Unknown encoding: {encoding}.")





