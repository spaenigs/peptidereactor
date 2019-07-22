from sklearn.externals import joblib as jl
import sys


sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


rule generate_lambda_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/apaac/csv/original/" + \
            "{dataset}_{part}_ifeature_apaacencoder_lambda-{lambdaValue}.csv"
    run:
        import encoder.ifeature.param_lambda.encoder as param_lambda_encoder
        for df in param_lambda_encoder.APAACEncoder(
                in_data=jl.load(str(input)),  cores=4,
                lambdaValue=int(wildcards.lambdaValue)).encode():
            df.to_csv(str(output))