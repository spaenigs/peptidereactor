from sklearn.manifold import TSNE

import pandas as pd
import numpy as np

import re


def run_clustering(pattern, smk_in, smk_out):

    df = pd.read_csv(smk_in, index_col=0)

    if not df.empty:

        tril_corr = df.applymap(lambda x: 1.0 if np.isnan(x) else x)

        res = tril_corr * tril_corr.transpose()

        embedding = TSNE(n_components=2, metric="precomputed")

        X_transformed = embedding.fit_transform(res)
        df_for_cls = pd.DataFrame(
            X_transformed,
            index=res.index,
            columns=["x1", "x2"])

        df_for_cls["type"] = list(map(lambda x: re.match(pattern, x).group(1), res.index))
        df_for_cls.to_csv(str(smk_out))

    else:
        df.to_csv(str(smk_out))