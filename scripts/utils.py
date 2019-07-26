from snakemake.io import expand
from encoder.ifeature.param_free.encoder import AAIndexEncoder

AAINDEX   = "aaindex"
AAC       = "aac"
APAAC     = "apaac"
BINARY    = "binary"
BLOSUM62  = "blosum62"
CKSAAGP   = "cksaagp"
CKSAAP    = "cksaap"  # thousands of columns and very sparse, no dataset passes rule filter_datasets
CTRIAD    = "ctriad"  # very sparse, no dataset passes rule filter_datasets
CTDT      = "ctdt"
CTDC      = "ctdc"
CTDD      = "ctdd"
DISORDER  = "disorder"
DDE       = "dde"
DPC       = "dpc"
EAAC      = "eaac"
EGAAC     = "egaac"
GAAC      = "gaac"
GEARY     = "geary"
GTPC      = "gtpc"
GDPC      = "gdpc"
KSCTRIAD  = "ksctriad"  # thousands of columns and very sparse, no dataset passes rule filter_datasets
MORAN     = "moran"
NMBROTO   = "nmbroto"
PAAC      = "paac"
PSEKRAAC  = "psekraac"
PSIPRED   = "psipred"
PSSM      = "pssm"
QSORDER   = "qsorder"
SOCNUMBER = "socnumber"
SPINEX    = "spinex"
TPC       = "tpc"
ZSCALE    = "zscale"

PARAM_BASED_ENCODINGS = ["apaac", "paac", "cksaagp", "cksaap", "ctriad",
                         "ksctriad", "geary", "moran", "nmbroto", "qsorder",
                         "socnumber", "eaac", "egaac"]

PARAM_FREE_ENCODINGS = ["binary", "aac", "gaac", "ctdt", "ctdc", "ctdd", "tpc",
                        "gtpc", "gtpc", "dpc", "gdpc", "dde", "blosum62", "zscale"]

REST_ENCODINGS = ["AAINDEX", "PSEKRAAC"]


ENCODING_PATTERN = {
    AAINDEX:   "aaindexencoder_aaindex-(.*?)\d+",
    APAAC:     "apaacencoder_(.*)",
    CKSAAGP:   "cksaagpencoder_(.*)",
    CKSAAP:    "cksaapencoder_(.*)",
    CTRIAD:    "ctriadencoder_(.*)",
    EAAC:      "eaacencoder_(.*)",
    EGAAC:     "egaacencoder_(.*)",
    GEARY:     "gearyencoder_(.*)",
    KSCTRIAD:  "ksctriadencoder_(.*)",
    MORAN:     "moranencoder_(.*)",
    NMBROTO:   "nmbrotoencoder_(.*)",
    PAAC:      "paacencoder_(.*)",
    QSORDER:   "qsorderencoder_(.*)",
    SOCNUMBER: "socnumberencoder_(.*)",
    PSEKRAAC:  "psekraac(.*?)_subtype"
}


def get_type(encoding, config):

    if encoding in [APAAC, PAAC]:
        return expand("{name}encoder_lambda-{lambdaValue}",
                      name=encoding,
                      lambdaValue=config["lambda_based"][encoding]["lambdas"])

    elif encoding in [CKSAAGP, CKSAAP, CTRIAD, KSCTRIAD]:
        return expand("{name}encoder_gap-{gapValue}",
                      name=encoding,
                      gapValue=config["gap_based"][encoding]["gaps"])

    elif encoding in [GEARY, MORAN, NMBROTO, QSORDER, SOCNUMBER]:
        return expand("{name}encoder_nlag-{nlag}",
               name=encoding,
               nlag=config["nlag_based"][encoding]["nlags"])

    elif encoding in [EAAC, EGAAC]:
        return expand("{name}encoder_window-{window}",
                      name=encoding,
                      window=config["window_based"][encoding]["windows"])

    elif encoding == PSEKRAAC:
        files = []
        for type_ in config["psekraac"]["types"]:
            files += expand("{name}_subtype-{subtype}_raactype-{raactype}_ktuple-{ktuple}_glValue-{glambda}",
                            name=config["psekraac"][type_]["name"],
                            subtype=config["psekraac"][type_]["subtypes"],
                            raactype=config["psekraac"][type_]["raactypes"],
                            ktuple=config["psekraac"][type_]["ktuples"],
                            glambda=config["psekraac"][type_]["glambdas"])
        return files

    elif encoding == AAINDEX:
        return expand("{name}encoder_aaindex-{aaindex}",
                      name=encoding,
                      aaindex=AAIndexEncoder._aaindex_names)

    elif encoding in PARAM_FREE_ENCODINGS:
        return f"{encoding}encoder"

    elif encoding == PSSM:
        return f"{encoding}encoder_pssm"

    elif encoding in [DISORDER, SPINEX, PSIPRED]:
        return expand("{name}encoder_{type_name}",
                      name=encoding,
                      type_name=config[f"{encoding}_based"]["names"])

    else:
        raise ValueError(f"Unknown encoding: {encoding}.")


def get_unique_types(encoding):

    if encoding in [APAAC, PAAC]:
        return ["lambda"]

    elif encoding in [CKSAAGP, CKSAAP, CTRIAD, KSCTRIAD]:
        return ["gap"]

    elif encoding in [EAAC, EGAAC]:
        return ["window"]

    elif encoding in [GEARY, MORAN, NMBROTO, QSORDER, SOCNUMBER]:
        return ["nlag"]

    elif encoding == PSEKRAAC:
        return ["type1", "type2", "type3A", "type3B", "type4", "type5" ,"type6A" ,"type6B",
                "type6C", "type7", "type8", "type9", "type10", "type11","type12", "type13",
                "type14", "type15", "type16"]

    elif encoding == AAINDEX:
        return ["aaindex"]

    else:
        raise ValueError(f"Unknown encoding: {encoding}.")


def determine_input(wildcards, config):

    if wildcards.encoding == AAINDEX:
        return "00_data/out/{dataset}/{dataset}_{part}/encodings/aaindex/" + \
               "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"

    elif wildcards.encoding in [APAAC, PAAC, CKSAAGP, CKSAAP, CTRIAD, KSCTRIAD,
                                GEARY, MORAN, NMBROTO, QSORDER, SOCNUMBER,
                                EAAC, EGAAC, PSEKRAAC]:
        return "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
               "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"

    elif wildcards.encoding in [BINARY, AAC, GAAC, CTDT, CTDC, CTDD, TPC, GTPC,
                                DPC, GDPC, DDE, BLOSUM62, ZSCALE, DISORDER,
                                SPINEX, PSSM, PSIPRED]:
        return expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
                      "{dataset}_{part}_{type}.csv",
                      dataset=wildcards.dataset,
                      part=wildcards.part,
                      encoding=wildcards.encoding,
                      type=get_type(wildcards.encoding, config))

    else:
        raise ValueError(f"Unknown encoding: {wildcards.encoding}.")
