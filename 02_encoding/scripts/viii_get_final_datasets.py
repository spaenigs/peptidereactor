from shutil import copyfile

import pandas as pd

import glob
import os
import sys

sys.path.append(os.getcwd())
sys.path.append(os.getcwd() + "/apps")

import scripts.utils as utils
from ix_get_aaindex_files import get_aaindex_files

def write_empty_file(path):
    with open(path, mode="w") as f:
        f.write("")
        f.flush()

def write_and_copy(sources_destination, output_path):
    with open(output_path, mode="a") as f:
        for src, dst in sources_destination:
            copyfile(src, dst)
            f.write(f"{src},{dst}\n")
            f.flush()

dataset = snakemake.wildcards.dataset
part = snakemake.wildcards.part
encoding = snakemake.wildcards.encoding
normalized = snakemake.wildcards.normalized

if encoding == utils.AAINDEX:

    df = pd.read_csv(str(snakemake.input), index_col = 0)
    res = get_aaindex_files(df)

    srcs = [glob.glob(f"00_data/out/{dataset}/{dataset}_{part}/encodings/" +
           f"aaindex/csv/normalized/{dataset}_{part}_{r}.csv")[0]
           for r in res]

    dsts = [f"00_data/out/{dataset}/{dataset}_{part}/encodings/aaindex/csv/final/" + \
            f"geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"
            for src in srcs]

    write_and_copy(list(zip(srcs, dsts)), str(snakemake.output))

elif encoding in [utils.APAAC, utils.PAAC, utils.CKSAAGP, utils.CKSAAP, utils.CTRIAD, utils.KSCTRIAD,
                  utils.GEARY, utils.MORAN, utils.NMBROTO, utils.QSORDER, utils.SOCNUMBER,
                  utils.EAAC, utils.EGAAC]:

    df_gm = pd.read_csv(str(snakemake.input), index_col=0)

    if not df_gm.empty:

        df_filtered = df_gm.loc[df_gm["min_gm"].apply(lambda x: x == True), :]
        src = glob.glob(
            f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
            f"*{df_filtered.name.values[0]}.csv")[0]
        dst = f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
              f"geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"
        write_and_copy(list(zip([src], [dst])), str(snakemake.output))

    else:
        write_empty_file(str(snakemake.output))

elif encoding == utils.PSEKRAAC:

    df_gm = pd.read_csv(str(snakemake.input), index_col=0)

    if not df_gm.empty:

        df_filtered = df_gm.loc[df_gm["min_gm"].apply(lambda x: x == True), :]
        src = [glob.glob(f"00_data/out/{dataset}/{dataset}_{part}/encodings/" + \
                         f"{encoding}/csv/original/*{t}.csv")[0]
               for t in df_filtered.name.values]
        dst = [f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
               f"geom_median/tsne/normalized-{normalized}/{os.path.basename(s)}"
               for s in src]
        write_and_copy(list(zip(src, dst)), str(snakemake.output))

    else:
        write_empty_file(str(snakemake.output))

elif encoding in [utils.BINARY, utils.AAC, utils.GAAC, utils.CTDT, utils.CTDC, utils.CTDD,
                  utils.TPC, utils.DPC, utils.GDPC, utils.DDE, utils.BLOSUM62, utils.ZSCALE,
                  utils.GTPC, utils.PSSM]:

    src = glob.glob(
        f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/*.csv")[0]
    dst = f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
          f"geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"
    write_and_copy(list(zip([src], [dst])), str(snakemake.output))

elif encoding in [utils.DISORDER, utils.SPINEX, utils.PSIPRED]:

    srcs = glob.glob(
        f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/*.csv")

    dsts = [f"00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" + \
            f"geom_median/tsne/normalized-{normalized}/{os.path.basename(src)}"
            for src in srcs]

    write_and_copy(list(zip(srcs, dsts)), str(snakemake.output))

else:
    raise ValueError(f"Unknown encoding: {encoding}.")




