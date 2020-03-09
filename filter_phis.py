import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool


def f(i):

    phi_dir = "data/bachem/machine_learning/phi_correlation/"

    df = pd.read_csv("data/bachem/machine_learning/top_encodings.csv", index_col=0) \
             .filter(regex="window_length_" + str(i), axis=0) \
             .sort_values(by="0", ascending=False) \
             .iloc[:100, :]

    file = f"phi_correlation_wl_{str(i)}.csv"
    df_phi = pd.read_csv(phi_dir + file, index_col=0)
    df_phi.index = pd.MultiIndex.from_arrays([df_phi["e1"], df_phi["e2"], df_phi["e3"]])
    df_phi["keep"] = False

    for i, series in df_phi.iterrows():
        if series["e1"] in df.index and \
                series["e2"] in df.index and \
                series["e3"] in df.index:
            df_phi.loc[(series["e1"], series["e2"], series["e3"]), "keep"] = True
            # for _i, _j, _k in permutations([series["e1"], series["e2"], series["e3"]]):
            #     if (_i, _j, _k) in df_phi.index:
            #         df_phi.loc[(_i, _j, _k), "keep"] = True

    # print(df_phi.shape)
    # print(df_phi.loc[df_phi["keep"] == True, :])

    df_phi.loc[df_phi["keep"] == True, :].to_csv(file)


p = Pool(6)
p.map(f, [8, 11, 13, 15, 17, 20])
