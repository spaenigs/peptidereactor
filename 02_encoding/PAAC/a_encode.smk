from sklearn.externals import joblib as jl
import sys


sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


localrules: paac_generate_lambda_based_encodings


rule paac_generate_lambda_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,paac}/csv/original/" + \
            "{dataset}_{part}_ifeature_{encoding}encoder_lambda-{lambdaValue}.csv"
    run:
        import encoder.ifeature.param_lambda.encoder as param_lambda_encoder
        for df in param_lambda_encoder.PAACEncoder(
                in_data=jl.load(str(input)),  cores=4,
                lambdaValue=int(wildcards.lambdaValue)).encode():
            df.to_csv(str(output))