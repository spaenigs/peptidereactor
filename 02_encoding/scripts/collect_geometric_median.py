import pandas as pd
from pandas.errors import EmptyDataError

dfs = list(snakemake.input)
res = pd.DataFrame()
for path in dfs:
    try:
        res = pd.concat([res, pd.read_csv(path, index_col=0)], axis=0)
    except EmptyDataError:
        pass
res.to_csv(str(snakemake.output))