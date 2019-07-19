from pandas.errors import EmptyDataError
from sklearn.preprocessing import MinMaxScaler

import pandas as pd

try:
    ds1 = pd.read_csv(str(snakemake.input), engine="c", index_col=0)
    if snakemake.config["normalize"] == "yes":
        vals = ds1.iloc[:, :-1].astype("float64").values
        ds1_normalized = pd.DataFrame(MinMaxScaler().fit_transform(vals),
                     index=ds1.index, columns=ds1.columns[:-1])
        ds1_normalized["y"] = ds1["y"]
        ds1_normalized.to_csv(str(snakemake.output))
    else:
        ds1.to_csv(str(snakemake.output))

except EmptyDataError as e:
    with open(str(snakemake.output), mode="w") as f:
        f.write("")
        f.flush()
