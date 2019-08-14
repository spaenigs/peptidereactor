import sys
import os
sys.path.append(os.getcwd())
sys.path.append(os.getcwd() + "/apps")

import scripts.utils as utils

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

res = pd.read_csv(str(snakemake.input), index_col=0)

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

plt.title(f"{utils.get_encoding_description(snakemake.wildcards.encoding)} ({snakemake.wildcards.encoding.upper()}):\n" +
          "paired t-test with corrected variance on classification error")

plt.savefig(str(snakemake.output), bbox_inches="tight")