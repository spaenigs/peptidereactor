from sklearn.externals import joblib as jl
import sys


sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


localrules: cksaagp_generate_encodings


rule cksaagp_generate_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,cksaagp}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_gap-{gapValue}.csv"
    run:
        import encoder.ifeature.param_gap.encoder as param_gap_encoder
        for df in param_gap_encoder.CKSAAGPEncoder(
                in_data=jl.load(str(input)),  cores=4,
                gap=int(wildcards.gapValue)).encode():
            df.to_csv(str(output))