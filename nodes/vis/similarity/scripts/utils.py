from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

import numpy as np


def path(paths, prop1, prop2): \
    return [p for p in paths if prop1 in p and prop2 in p][0]


def cluster(values, axis):
    if axis == 1:
        values = values.T
    linkage = hierarchy.linkage(pdist(values), method="average", metric="euclidean")
    return hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)["leaves"]
