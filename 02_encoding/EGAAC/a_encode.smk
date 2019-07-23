from sklearn.externals import joblib as jl
import sys


sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


localrules: egaac_generate_encodings


rule egaac_generate_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,egaac}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_window-{windowValue}.csv"
    run:
        import encoder.ifeature.param_window.encoder as param_window_encoder
        for df in param_window_encoder.EGAACEncoder(
                in_data=jl.load(str(input)),  cores=4,
                window=int(wildcards.windowValue)).encode():
            df.to_csv(str(output))