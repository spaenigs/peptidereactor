from sklearn.externals import joblib as jl

sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])


localrules: generate_psekraac_based_encodings


rule generate_psekraac_based_encodings:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/psekraac/csv/original/" + \
            "{dataset}_{part}_ifeature_{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}.csv"
    run:
        import encoder.ifeature.psekraac.encoder as psekraac_encoder
        if wildcards.name == "psekraactype1encoder":
            for df in psekraac_encoder.PseKRAACType1Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype2encoder":
            for df in psekraac_encoder.PseKRAACType2Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype3Aencoder":
            for df in psekraac_encoder.PseKRAACType3AEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype3Bencoder":
            for df in psekraac_encoder.PseKRAACType3BEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype4encoder":
            for df in psekraac_encoder.PseKRAACType4Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype5encoder":
            for df in psekraac_encoder.PseKRAACType5Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype6Aencoder":
            for df in psekraac_encoder.PseKRAACType6AEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype6Bencoder":
            for df in psekraac_encoder.PseKRAACType6BEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype6Cencoder":
            for df in psekraac_encoder.PseKRAACType6CEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype7encoder":
            for df in psekraac_encoder.PseKRAACType7Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype8encoder":
            for df in psekraac_encoder.PseKRAACType8Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype9encoder":
            for df in psekraac_encoder.PseKRAACType9Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype10encoder":
            for df in psekraac_encoder.PseKRAACType10Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype11encoder":
            for df in psekraac_encoder.PseKRAACType11Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype12encoder":
            for df in psekraac_encoder.PseKRAACType12Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype13encoder":
            for df in psekraac_encoder.PseKRAACType13Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype14encoder":
            for df in psekraac_encoder.PseKRAACType14Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype15encoder":
            for df in psekraac_encoder.PseKRAACType15Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        else:
            for df in psekraac_encoder.PseKRAACType16Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))

