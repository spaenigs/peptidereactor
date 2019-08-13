import scripts.utils as utils
import pandas as pd
import numpy as np
import itertools


rule t_test_on_scores_error:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
         "{dataset}_{part}_normalized-{normalized}_cross_validation_scores.csv"
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
          "{dataset}_{part}_normalized-{normalized}_ttest_error_scores.csv"
    run:
        from scipy import stats
        df = pd.read_csv(str(input), index_col=0)
        names = df.index
        res = pd.DataFrame(np.zeros((len(names), len(names))) + 1.0,
                           index=names,
                           columns=names)
        for n1, n2 in itertools.combinations(names, 2):
            _, pvalue = stats.ttest_rel(df.loc[n1, :].values, df.loc[n1, :].values)
            res.loc[n1, n2] = pvalue
        res = res.transpose() * res
        res.to_csv(str(output))


rule plot_t_test_on_scores_error:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/cv/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error_scores.csv"
    output:
        "00_data/out/neuropeptides/plots/{encoding}/" + \
        "{dataset}_{part}_normalized-{normalized}_ttest_error_scores.pdf"
    run:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        res = pd.read_csv(str(input), index_col=0)

        # https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
        fig, ax = plt.subplots()
        pc = ax.imshow(res.values, vmin=0, vmax=1, cmap='Greys')
        ax.set_xticks(np.arange(len(res.columns)))
        ax.set_yticks(np.arange(len(res.index)))
        ax.set_xticklabels(res.columns)
        ax.set_yticklabels(res.index)

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        cbar = ax.figure.colorbar(pc)
        cbar.ax.set_ylabel("p-value, (*): <= 0.05, (**): <= 0.01", rotation=-90, va="bottom")

        for i in range(len(res.index)):
            for j in range(len(res.columns)):
                if res.loc[res.index[i], res.columns[j]] <= 0.01:
                    ax.text(j, i, "(**)", ha="center", va="center", color="black")
                elif res.loc[res.index[i], res.columns[j]] <= 0.05:
                    ax.text(j, i, "(*)", ha="center", va="center", color="black")

        fig.tight_layout()
        fig.set_size_inches(9, 7)

        plt.title(f"{utils.get_encoding_description(wildcards.encoding)} ({wildcards.encoding.upper()}):\n" +
                  "paired t-test with corrected variance on classification error")

        plt.savefig(str(output), bbox_inches="tight")