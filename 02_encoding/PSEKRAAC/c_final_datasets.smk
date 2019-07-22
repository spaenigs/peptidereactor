import itertools

from pandas.errors import EmptyDataError
from typing import List

import numpy as np
import pandas as pd
import glob, re, os

from snakemake.io import expand, temp


def write_empty_file(path):
    with open(path, mode="w") as f:
        f.write("")
        f.flush()


rule generate_distance_matrix:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/csv/normalized/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}_normalized-{normalized}.csv",
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/csv/normalized/{dataset}_{part}_normalized-{normalized}.txt"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/correlation/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}_normalized-{normalized}_vs_rest.csv"
    run:
        from scipy import interpolate, stats

        def interpolate_to(X: np.array, dim: int):
            ydim = X.shape[1]
            x = np.arange(0, ydim if ydim >= 2 else ydim + 1)
            func = interpolate.interp1d(x, X, kind="linear")
            return func(np.linspace(0, ydim-1, dim))

        try:
            actual_csv = str(input[0])
            ds1 = pd.read_csv(actual_csv, engine="c", index_col=0).iloc[:, :-1].astype("float64").values

            files = [f for f in
                     glob.glob(f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/encodings/psekraac/csv/normalized/" + \
                                    f"{wildcards.dataset}_{wildcards.part}_ifeature_*_normalized-{wildcards.normalized}*")
                     # glob.glob(f"data/out/{wildcards.encoding}/csv/normalized/{wildcards.dataset}_{wildcards.source}_*.normalized-{wildcards.normalized}*")
                     if os.path.getsize(f) > 0]

            files_filtered =\
                list(map(lambda tup: tup[1], filter(lambda tup: tup[0] == actual_csv, itertools.combinations(files, 2))))

            # pattern = f"data/out/{wildcards.encoding}/csv/normalized/{wildcards.dataset}_{wildcards.source}_(.*?)\.normalized-{wildcards.normalized}\.csv"
            pattern = f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/encodings/psekraac/csv/normalized/" + \
                            f"{wildcards.dataset}_{wildcards.part}_ifeature_(.*?)_normalized-{wildcards.normalized}.csv"

            row_name = re.match(pattern, actual_csv).group(1)

            inner_dict = {}
            for _f in files:
                inner_dict.setdefault(re.match(pattern, _f).group(1), "NA")

            res_corr = {row_name: inner_dict}

            if len(files_filtered) is not 0:
                for ff in files_filtered:
                    col_name = re.match(pattern, ff).group(1)
                    ds2 = pd.read_csv(ff, engine="c", index_col=0).iloc[:, :-1].values
                    to_dim = int(np.ceil(np.mean([ds1.shape[1], ds2.shape[1]])))
                    ds1_interpolated, ds2_interpolated = interpolate_to(ds1, to_dim), interpolate_to(ds2, to_dim)
                    corr, _ = stats.pearsonr(ds1_interpolated.ravel(), ds2_interpolated.ravel())
                    res_corr[row_name][col_name] = 1 - corr

            df = pd.DataFrame(res_corr)
            df.to_csv(str(output))

        except EmptyDataError as e:
            write_empty_file(str(output))


def collect_files(wildcards):
    files = []
    for type_ in config["psekraac"]["types"]:
        files += expand("00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/correlation/" + \
                            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}_normalized-{normalized}_vs_rest.csv",
                        dataset=wildcards.dataset, part=wildcards.part,  normalized=wildcards.normalized,
                        name=config["psekraac"][type_]["name"],
                        subtype=config["psekraac"][type_]["subtypes"],
                        raactype=config["psekraac"][type_]["raactypes"],
                        ktuple=config["psekraac"][type_]["ktuples"],
                        glambda=config["psekraac"][type_]["glambdas"])
    return files

rule collect_distance_matrix:
    input:
        collect_files
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    run:
         def sort_df_to_lower_tril(df):
            """
            1) Sort rows
                a   b    c              a   b  c
            a  NaN NaN  1.0         b  NaN NaN  NaN
            b  NaN NaN  NaN  ==>    a  NaN NaN  1.0
            c  3.0 NaN  2.0         c  3.0 NaN  2.0
            2) Sort cols
                a   b  c                c   a   b
            b  NaN NaN  NaN         b  NaN NaN NaN
            a  NaN NaN  1.0  ==>    a  1.0 NaN NaN
            c  3.0 NaN  2.0         c  2.0 3.0 NaN

            """
            sorted_row_names = np.array(sorted([[row[0], ~np.isnan(row[1]).sum()]
                                                for row in df.iterrows()],
                                               key=lambda tup: tup[1]))[:,0]
            df = df.reindex(sorted_row_names)

            sorted_col_names = np.array(sorted([[row[0], np.isnan(row[1]).sum()]
                                                for row in df.transpose().iterrows()],
                                               key=lambda tup: tup[1]))[:,0]
            return df[sorted_col_names]

         def concat(paths: List[str]) -> pd.DataFrame:
             df_res = pd.DataFrame()
             for path in paths:
                 try:
                    df = pd.read_csv(path, index_col=0)
                    df_res = pd.concat([df_res, df], axis=1, sort=True)
                 except EmptyDataError:
                     pass
             return df_res if df_res.empty else sort_df_to_lower_tril(df_res)

         concat(list(input)).to_csv(str(output))


rule run_clustering:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized,yes|no}.csv"
    run:
        from sklearn.manifold import MDS, TSNE

        df = pd.read_csv(str(input), index_col=0)

        if not df.empty:

            tril_corr = df.applymap(lambda x: 1.0 if np.isnan(x) else x)

            res = tril_corr * tril_corr.transpose()

            embedding = TSNE(n_components=2, metric="precomputed")

            X_transformed = embedding.fit_transform(res)
            df_for_cls = pd.DataFrame(
                X_transformed,
                index=res.index,
                columns=["x1", "x2"])

            df_for_cls["type"] = list(map(lambda x: re.match("psekraac(.*?)_subtype", x).group(1), res.index))
            df_for_cls.to_csv(str(output))

        else:
            df.to_csv(str(output))


rule compute_geometric_median:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/geom_median/{dataset}_{part}_normalized-{normalized}_{type}_vs_rest.csv"
    run:
        from scipy.spatial.distance import euclidean

        df = pd.read_csv(str(input), index_col=0)
        df["name"] = df.index

        df_type_filtered = df.loc[df.name.apply(lambda x: wildcards.type in x), :]

        names = df_type_filtered.name
        points = dict(zip(names, map(lambda x: dict(zip(names, map(lambda x: 0.0, range(len(names))))),
                                     range(len(names)))))

        res = {}
        for i in names:
            for j in names:
                if i != j:
                    if points[i][j] == 0.0:
                        # https://en.wikipedia.org/wiki/Geometric_median
                        points[i][j] = euclidean(df_type_filtered.loc[i, ["x1", "x2"]], df_type_filtered.loc[j, ["x1", "x2"]])
                    else:
                        continue
            res[i] = [sum(sorted(points[i].values()))]

        res_df = pd.DataFrame(res).transpose()

        if not res_df.empty:
            res_df.columns = ["distance_gm"]

            tmp2 = pd.concat([df_type_filtered, res_df], axis=1)
            tmp2["min_gm"] = False

            actual_total_min = 1e6
            actual_total_min_name = ""
            for row in tmp2.loc[tmp2.name.apply(lambda x: wildcards.type in x), :].iterrows():
                if row[1]["distance_gm"] < actual_total_min:
                   actual_total_min = row[1]["distance_gm"]
                   actual_total_min_name = row[0]
            tmp2.loc[actual_total_min_name, "min_gm"] = True

            tmp2.to_csv(str(output))

        else:
            write_empty_file(str(output))


rule collect_geometric_median:
    input:
        lambda wildcards: expand("00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/geom_median/" + \
                                    "{dataset}_{part}_normalized-{normalized}_{type}_vs_rest.csv",
                                 dataset=wildcards.dataset, part=wildcards.part,
                                 normalized=wildcards.normalized,
                                 type=[f"{t}encoder" for t in config["psekraac"]["types"]])
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    run:
        dfs = list(input)
        res = pd.DataFrame()
        for path in dfs:
            try:
                res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0)
            except EmptyDataError:
                pass
        res.to_csv(str(output))


rule plot_clustering:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        "00_data/out/{dataset}/plots/{dataset}_{part}_psekraac_normalized-{normalized}_tsne.svg"
    run:
        import matplotlib.pyplot as plt
        from matplotlib import colors as mcolors

        colors = list(dict(mcolors.TABLEAU_COLORS, **mcolors.CSS4_COLORS).keys())

        fig, ax = plt.subplots()

        df_gm = pd.read_csv(str(input), index_col=0)

        if not df_gm.empty:

            unames = df_gm.type.unique()

            for i in range(len(unames)):
                df_filtered = df_gm.loc[df_gm.name.apply(lambda x: unames[i] in x), :]
                ax.scatter(df_filtered.x1, df_filtered.x2, c=colors[i], label=f"{unames[i]}", alpha=0.2, marker=".")

            for i in range(len(unames)):
                df_filtered = df_gm.loc[df_gm.name.apply(lambda x: unames[i] in x), :]
                df_filtered = df_filtered.loc[df_filtered["min_gm"].apply(lambda x: x == True), :]
                ax.scatter(df_filtered.x1, df_filtered.x2, marker="s", s=70, label=f"(geom. median)", c=colors[i])
                if not df_filtered.empty:
                    ax.text(df_filtered.x1,
                            df_filtered.x2 - (0.2 * df_filtered.x2),
                            df_filtered.name.values[0],
                            ha="center")

            ax.legend(ncol=3, bbox_to_anchor=(1, 1), loc="upper left")
            ax.grid(True)

            fig = plt.gcf()
            fig.set_size_inches(20, 13)

            plt.subplots_adjust(right=0.75)
            plt.title(f"{wildcards.dataset}{', normalized, ' if wildcards.normalized == 'yes' else ', ' }{wildcards.clustering.upper()}")
            plt.savefig(str(output))

        else:
            plt.savefig(str(output))


rule get_final_datasets:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        temp("00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/tsne/{dataset}_{part}_normalized-{normalized}_final_datasets.txt")
    run:
        from shutil import copyfile

        df_gm = pd.read_csv(str(input), index_col=0)

        if (not df_gm.empty) or (not df_gm.empty):

            with open(str(output), mode="a") as f:

                df_filtered = df_gm.loc[df_gm["min_gm"].apply(lambda x: x == True), :]
                src = glob.glob(f"data/in/csv/{wildcards.encoding}/*{df_filtered.name.values[0]}.csv")[0]
                dst = f"data/out/{wildcards.encoding}/csv/final/{wildcards.kind}/{wildcards.clustering}/normalized-{wildcards.normalized}/{os.path.basename(src)}"
                copyfile(src, dst)
                f.write(os.path.basename(src))

                f.flush()

        else:
            write_empty_file(str(output))