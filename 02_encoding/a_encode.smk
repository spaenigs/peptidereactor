from sklearn.externals import joblib as jl

import os
os.environ["R_LIBS_USER"] = config["cwd"] + "/" + config["programs"]["rlibs"]

import sys
sys.path.append(config["cwd"] + "/" + config["programs"]["iFeature"])

try:
    import encoder.ifeature.param_free.encoder as param_free_encoder
except ModuleNotFoundError as e:
    print()
    print(e)
    print("Loaded environment specific configuration?", file=sys.stderr)
    sys.exit(1)


localrules: generate_psekraac_based_encodings


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


rule generate_binary_encoding:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/original/" + \
        "{dataset}_{part}_binaryencoder"
    run:
        df = param_free_encoder.BinaryEncoder(
            in_data=jl.load(str(input)), cores=4, run_msa=False
        ).encode()
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
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib"
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
        else:
            df = param_free_encoder.ZscaleEncoder(in_data=jl.load(str(input)), cores=4)\
                .encode()
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




### TODO ####


# rule annotate_sequence_names:
#     input:
#         "data/out/joblib/{dataset}-pssms_filtered.joblib",
#         "data/out/joblib/{dataset}-pssms_filtered-msa.joblib"
#     output:
#         "data/out/joblib/{dataset,[A-Za-z]+(_ds[12])?}-annotated.joblib",
#         "data/out/joblib/{dataset,[A-Za-z]+(_ds[12])?}-annotated-msa.joblib"
#     run:
#         def add_names(input_data_):
#             res_seqs, res_classes = [], []
#             for tup in zip(*input_data_):
#                 seq_tup = tup[0]
# 		print(f"original seq: {seq_tup[0]}")
#                 print(f"with wildcards: {wildcards.dataset}-{seq_tup[0]}")
#                 seq_tup[0] = f"{wildcards.dataset}-{seq_tup[0]}"
#                 res_seqs.append(seq_tup)
#                 res_classes.append(tup[1])
#             return res_seqs, res_classes
#         input_data, input_data_msa = jl.load(str(input[0])), jl.load(str(input[1]))
#         jl.dump(value=add_names(input_data), filename=str(output[0]))
#         jl.dump(value=add_names(input_data_msa), filename=str(output[1]))
#
#
# def generate_vsl2_based_file_names(wildcards):
#     seq_names, files = get_sequence_names(wildcards), []
#     files += expand("data/out/profile/{dataset}-{seq_name}.dis",
#                     dataset=wildcards.dataset,
#                     seq_name=seq_names)
#     files += expand("data/out/profile/{dataset}-{seq_name}.flat",
#                     dataset=wildcards.dataset,
#                     seq_name=seq_names)
#     return files
#
#
# rule generate_disorder_based_encodings:
#     input:
#          "data/out/joblib/{dataset}-annotated.joblib",
#          "data/out/joblib/{dataset}-annotated-msa.joblib",
#          generate_vsl2_based_file_names
#     output:
#         "data/out/csv/disorder/{dataset}_ifeature_disorder_{name}-encoder.csv"
#     run:
#         import encoder.ifeature.pssm.disprot_encoding.encoder as disprot_encoder
#         if wildcards.encoding == "disorder":
#             df = disprot_encoder.DisorderEncoder(
#                 in_data=jl.load(str(input[0])),
#                 profile_dir="data/out/profile",
#                 cores=1
#             ).encode()
#         elif wildcards.encoding == "disorderb":
#             df = disprot_encoder.DisorderBEncoder(
#                 in_data=jl.load(str(input[1])),
#                 profile_dir="data/out/profile",
#                 cores=1, run_msa=False
#             ).encode()
#         else:
#             df = disprot_encoder.DisorderCEncoder(
#                 in_data=jl.load(str(input[0])),
#                 profile_dir="data/out/profile",
#                 cores=1
#             ).encode()
#         df.to_csv(str(output))
#
#
# def generate_spinex_based_file_names(wildcards):
#     seq_names, files = get_sequence_names(wildcards), []
#     files += expand("data/out/profile/{dataset}-{seq_name}.spXout",
#                     dataset=wildcards.dataset,
#                     seq_name=seq_names)
#     return files
#
#
# rule generate_spinex_based_encodings:
#     input:
#         "data/out/joblib/{dataset}-annotated.joblib",
#         generate_spinex_based_file_names
#     output:
#         "data/out/csv/spinex/{dataset}_ifeature_spinex_{name}-encoder.csv"
#     run:
#         import encoder.ifeature.pssm.spinex_encoding.encoder as spx_encoder
#         if wildcards.encoding == "asa":
#             df = spx_encoder.ASAEncoder(
#                 in_data=jl.load(str(input[0])),
#                 profile_dir="data/out/profile",
#                 cores=1,
#             ).encode()
#         else:
#             df = spx_encoder.TAEncoder(
#                 in_data=jl.load(str(input[0])),
#                 profile_dir="data/out/profile",
#                 cores=1,
#             ).encode()
#         df.to_csv(str(output))
#
#
# rule generate_pssm_encoding:
#     input:
#         "data/out/joblib/{dataset}-annotated.joblib"
#     output:
#         "data/out/csv/pssm/{dataset}_ifeature_pssm_{name}-encoder.csv"
#     run:
#         import encoder.ifeature.pssm.pssm_encoding.encoder as pssm_encoder
#         df = pssm_encoder.PSSMEncoder(
#             in_data=jl.load(str(input)),
#             profile_dir="data/out/profile",
#             cores=1
#         ).encode()
#         df.to_csv(str(output))
#
#
# def generate_psipred_based_fime_names(wildcards):
#     seq_names, files = get_sequence_names(wildcards), []
#     files += expand("data/out/profile/{dataset}-{seq_name}.spXout",
#                     dataset=wildcards.dataset,
#                     seq_name=seq_names)
#     files += expand("data/out/profile/{dataset}-{seq_name}.ss2",
#                     dataset=wildcards.dataset,
#                     seq_name=seq_names)
#     return files
#
#
# rule generate_psipred_based_encodings:
#     input:
#         "data/out/joblib/{dataset}-annotated.joblib",
#         "data/out/joblib/{dataset}-annotated-msa.joblib",
#         generate_psipred_based_fime_names
#     output:
#         "data/out/csv/psipred/{dataset}_ifeature_psipred_{name}-encoder.csv"
#     run:
#         import encoder.ifeature.pssm.psipred_encoding.encoder as psipred_encoder
#         if wildcards.encoding == "sseb":
#             df = psipred_encoder.SSEBEncoder(
#                 in_data=jl.load(str(input[1])),
#                 profile_dir="data/out/profile",
#                 cores=1, run_msa=False
#             ).encode()
#         else:
#             df = psipred_encoder.SSECEncoder(
#                 in_data=jl.load(str(input[0])),
#                 profile_dir="data/out/profile",
#                 cores=1,
#             ).encode()
#         df.to_csv(str(output))


