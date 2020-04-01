import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool


def f(i):

    phi_dir = "data/pds/machine_learning/phi_correlation/"

    df = pd.read_csv("data/pds/machine_learning/top_encodings.csv", index_col=0) \
             .filter(regex="window_length_" + str(i), axis=0) \
             .sort_values(by="0", ascending=False) \
             .iloc[:40, :]

    file = f"phi_correlation_wl_{str(i)}.csv"
    df_phi = pd.read_csv(phi_dir + file, index_col=0)
    df_phi.index = pd.MultiIndex.from_arrays([df_phi["e1"], df_phi["e2"], df_phi["e3"]])
    df_phi["keep"] = False

    for i, series in df_phi.iterrows():
        if "psekraac" in f"{series['e1']}{series['e2']}{series['e3']}":
            continue
        if series["e1"] in df.index and \
                series["e2"] in df.index and \
                series["e3"] in df.index:
            df_phi.loc[(series["e1"], series["e2"], series["e3"]), "keep"] = True

    df_phi.index = pd.Index(range(df_phi.shape[0]))
    df_phi.loc[df_phi["keep"] == True, :].to_csv(file)


p = Pool(6)
p.map(f, [8, 11, 13, 15, 17, 20])
