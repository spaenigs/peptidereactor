from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

import numpy as np


def is_struc_based(e):
    if "asa" in e:
        return True
    elif "delaunay" in e:
        return True
    elif "delaun" in e:
        return True
    elif "disorder" in e:
        return True
    elif "disord" in e:
        return True
    elif "elect" in e:
        return True
    elif "hull" in e:
        return True
    elif "qsar" in e:
        return True
    elif "sse" in e:
        return True
    elif "ta" in e:
        return True
    else:
        return False


def cluster(values, axis):
    if axis == 1:
        values = values.T
    linkage = hierarchy.linkage(pdist(values), method="average", metric="euclidean")
    return hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)["leaves"]
