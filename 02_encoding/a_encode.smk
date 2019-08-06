from sklearn.externals import joblib as jl

import os
os.environ["R_LIBS_USER"] = config["cwd"] + "/" + config["programs"]["rlibs"]

localrules: generate_window_based_encodings,
            generate_lambda_based_encodings,
            generate_gap_based_encodings,
            generate_nlag_based_encodings,
            generate_psekraac_based_encodings,
            generate_paramfree_based_encoding,
            generate_aaindex_encoding,
            generate_disorder_based_encodings,
            generate_psipred_based_encodings,
            generate_pssm_encoding

rule generate_window_based_encodings:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_{encoding}encoder_window-{windowValue}.csv"
    run:
        import encoder.ifeature.param_window.encoder as param_window_encoder
        if wildcards.encoding == "eaac":
            for df in param_window_encoder.EAACEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    window=int(wildcards.windowValue)).encode():
                df.to_csv(str(output))
        else:
            for df in param_window_encoder.EGAACEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    window=int(wildcards.windowValue)).encode():
                df.to_csv(str(output))


rule generate_lambda_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_{encoding}encoder_lambda-{lambdaValue}.csv"
    run:
        import encoder.ifeature.param_lambda.encoder as param_lambda_encoder
        if wildcards.encoding == "paac":
            for df in param_lambda_encoder.PAACEncoder(
                    in_data=jl.load(str(input)),  cores=4,
                    lambdaValue=int(wildcards.lambdaValue)).encode():
                df.to_csv(str(output))
        else:
            for df in param_lambda_encoder.APAACEncoder(
                    in_data=jl.load(str(input)),  cores=4,
                    lambdaValue=int(wildcards.lambdaValue)).encode():
                df.to_csv(str(output))


rule generate_nlag_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_{encoding}encoder_nlag-{nlagValue}.csv"
    run:
        import encoder.ifeature.param_nlag.encoder as param_nlag_encoder
        if wildcards.encoding == "moran":
            for df in param_nlag_encoder.MoranEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    nlag=int(wildcards.nlagValue)).encode():
                df.to_csv(str(output))
        elif wildcards.encoding == "geary":
            for df in param_nlag_encoder.GearyEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    nlag=int(wildcards.nlagValue)).encode():
                df.to_csv(str(output))
        elif wildcards.encoding == "nmbroto":
            for df in param_nlag_encoder.NMBrotoEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    nlag=int(wildcards.nlagValue)).encode():
                df.to_csv(str(output))
        elif wildcards.encoding == "socnumber":
            for df in param_nlag_encoder.SOCNumberEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    nlag=int(wildcards.nlagValue)).encode():
                df.to_csv(str(output))
        else:
            for df in param_nlag_encoder.QSOrderEncoder(
                    in_data=jl.load(str(input)), cores=4,
                    nlag=int(wildcards.nlagValue)).encode():
                df.to_csv(str(output))


rule generate_gap_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_{encoding}encoder_gap-{gapValue}.csv"
    run:
        import encoder.ifeature.param_gap.encoder as param_gap_encoder
        if wildcards.encoding == "cksaap":
                for df in param_gap_encoder.CKSAAPEncoder(
                        in_data=jl.load(str(input)), cores=4,
                        gap=int(wildcards.gapValue)).encode():
                    df.to_csv(str(output))
        elif wildcards.encoding == "cksaagpencoder":
                for df in param_gap_encoder.CKSAAGPEncoder(
                        in_data=jl.load(str(input)), cores=4,
                        gap=int(wildcards.gapValue)).encode():
                    df.to_csv(str(output))
        elif wildcards.encoding == "ctriadencoder":
                for df in param_gap_encoder.CTriadEncoder(
                        in_data=jl.load(str(input)), cores=4,
                        gap=int(wildcards.gapValue)).encode():
                    df.to_csv(str(output))
        else:
             for df in param_gap_encoder.KSCTriadEncoder(
                        in_data=jl.load(str(input)), cores=4,
                        gap=int(wildcards.gapValue)).encode():
                    df.to_csv(str(output))


rule generate_aaindex_encoding:
     input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
     output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_aaindexencoder_aaindex-{aaindex}.csv"
     run:
        for df in param_free_encoder.AAIndexEncoder(
                in_data=jl.load(str(input)), cores=4, index=wildcards.aaindex).encode():
            df.to_csv(str(output))


rule generate_paramfree_based_encoding:
    input:
        lambda wildcards: \
            f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/joblib/" + \
            f"{wildcards.dataset}_{wildcards.part}_pssms_filtered{'_msa' if wildcards.encoding == 'binary' else ''}.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_{name}encoder.csv"
    run:
        if wildcards.encoding == "gaac":
            df = param_free_encoder.GAACEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "aac":
            df = param_free_encoder.AACEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "ctdt":
            df = param_free_encoder.CTDTEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "ctdc":
            df = param_free_encoder.CTDCEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "ctdd":
            df = param_free_encoder.CTDDEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "tpc":
            df = param_free_encoder.TPCEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "gtpc":
            df = param_free_encoder.GTPCEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "dpc":
            df = param_free_encoder.DPCEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "gdpc":
            df = param_free_encoder.GDPCEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "dde":
            df = param_free_encoder.DDEEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "blosum62":
            df = param_free_encoder.Blosum62Encoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        elif wildcards.encoding == "binary":
            df = param_free_encoder.BinaryEncoder(
                in_data=jl.load(str(input)), cores=4, run_msa=False
            ).encode()
            df.to_csv(str(output))
        elif wildcards.encoding == "zscale":
            df = param_free_encoder.ZscaleEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
        else:
            raise ValueError("Unknown encoding.")
        df.to_csv(str(output))


rule generate_psekraac_based_encodings:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding,psekraac}/csv/original/" + \
         "{dataset}_{part}_{name}encoder_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}.csv"
    run:
        import encoder.ifeature.psekraac.encoder as psekraac_encoder
        if wildcards.name == "psekraactype1":
            for df in psekraac_encoder.PseKRAACType1Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype2":
            for df in psekraac_encoder.PseKRAACType2Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype3A":
            for df in psekraac_encoder.PseKRAACType3AEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype3B":
            for df in psekraac_encoder.PseKRAACType3BEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype4":
            for df in psekraac_encoder.PseKRAACType4Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype5":
            for df in psekraac_encoder.PseKRAACType5Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype6A":
            for df in psekraac_encoder.PseKRAACType6AEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype6B":
            for df in psekraac_encoder.PseKRAACType6BEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype6C":
            for df in psekraac_encoder.PseKRAACType6CEncoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype7":
            for df in psekraac_encoder.PseKRAACType7Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype8":
            for df in psekraac_encoder.PseKRAACType8Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype9":
            for df in psekraac_encoder.PseKRAACType9Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype10":
            for df in psekraac_encoder.PseKRAACType10Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype11":
            for df in psekraac_encoder.PseKRAACType11Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype12":
            for df in psekraac_encoder.PseKRAACType12Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype13":
            for df in psekraac_encoder.PseKRAACType13Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype14":
            for df in psekraac_encoder.PseKRAACType14Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype15":
            for df in psekraac_encoder.PseKRAACType15Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        elif wildcards.name == "psekraactype16":
            for df in psekraac_encoder.PseKRAACType16Encoder(
                in_data=jl.load(str(input)), cores=4,
                subtype=wildcards.subtype,
                raactype=int(wildcards.raactype),
                ktuple=int(wildcards.ktuple),
                glambda=int(wildcards.glambda)).encode():
                df.to_csv(str(output))
        else:
            raise ValueError(f"Unknown encoding: {wildcards.name}")


rule generate_disorder_based_encodings:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated.joblib",
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated_msa.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_disorderencoder_{name}.csv"
    run:
        import encoder.ifeature.pssm.disprot_encoding.encoder as disprot_encoder
        if wildcards.name == "disorder":
            df = disprot_encoder.DisorderEncoder(
                in_data=jl.load(str(input[0])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1
            ).encode()
        elif wildcards.name == "disorderb":
            df = disprot_encoder.DisorderBEncoder(
                in_data=jl.load(str(input[1])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1, run_msa=False
            ).encode()
        else:
            df = disprot_encoder.DisorderCEncoder(
                in_data=jl.load(str(input[0])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1
            ).encode()
        df.to_csv(str(output))


rule generate_spinex_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated.joblib",
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_spinexencoder_{name}.csv"
    run:
        import encoder.ifeature.pssm.spinex_encoding.encoder as spx_encoder
        if wildcards.name == "asa":
            df = spx_encoder.ASAEncoder(
                in_data=jl.load(str(input[0])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1,
            ).encode()
        else:
            df = spx_encoder.TAEncoder(
                in_data=jl.load(str(input[0])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1,
            ).encode()
        df.to_csv(str(output))


rule generate_pssm_encoding:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated.joblib",
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_pssmencoder_pssm.csv"
    run:
        import encoder.ifeature.pssm.pssm_encoding.encoder as pssm_encoder
        df = pssm_encoder.PSSMEncoder(
            in_data=jl.load(str(input)),
            profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
            cores=1
        ).encode()
        df.to_csv(str(output))


rule generate_psipred_based_encodings:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated.joblib",
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated_msa.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_psipredencoder_{name}.csv"
    run:
        import encoder.ifeature.pssm.psipred_encoding.encoder as psipred_encoder
        if wildcards.encoding == "sseb":
            df = psipred_encoder.SSEBEncoder(
                in_data=jl.load(str(input[1])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1, run_msa=False
            ).encode()
        else:
            df = psipred_encoder.SSECEncoder(
                in_data=jl.load(str(input[0])),
                profile_dir=f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
                cores=1,
            ).encode()
        df.to_csv(str(output))


