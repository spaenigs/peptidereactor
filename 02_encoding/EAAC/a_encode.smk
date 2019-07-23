from sklearn.externals import joblib as jl
import sys


sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


localrules: eaac_generate_encodings


rule eaac_generate_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,eaac}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_window-{windowValue}.csv"
    run:
        import encoder.ifeature.param_window.encoder as param_window_encoder
        for df in param_window_encoder.EAACEncoder(
                in_data=jl.load(str(input)),  cores=4,
                window=int(wildcards.windowValue)).encode():
            df.to_csv(str(output))