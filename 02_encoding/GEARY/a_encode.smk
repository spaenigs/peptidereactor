from sklearn.externals import joblib as jl
import sys


sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


localrules: geary_generate_encodings


rule geary_generate_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,geary}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_nlag-{nlagValue}.csv"
    run:
        import encoder.ifeature.param_nlag.encoder as param_nlag_encoder
        for df in param_nlag_encoder.GearyEncoder(
                in_data=jl.load(str(input)),  cores=4,
                nlag=int(wildcards.nlagValue)).encode():
            df.to_csv(str(output))