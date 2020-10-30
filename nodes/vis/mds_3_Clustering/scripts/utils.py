from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

import numpy as np


def cluster(values, axis):
    if axis == 1:
        values = values.T
    linkage = hierarchy.linkage(pdist(values), method="average", metric="euclidean")
    return hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)["leaves"]
