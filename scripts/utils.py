from snakemake.io import expand
from encoder.ifeature.param_free.encoder import AAIndexEncoder

AAC       = "aac"
AAINDEX   = "aaindex"
APAAC     = "apaac"
BINARY    = "binary"
BLOSUM62  = "blosum62"
CKSAAGP   = "cksaagp"
CKSAAP    = "cksaap"  # thousands of columns and very sparse, no dataset passes rule filter_datasets
CTDC      = "ctdc"
CTDD      = "ctdd"
CTDT      = "ctdt"
CTRIAD    = "ctriad"  # very sparse, no dataset passes rule filter_datasets
DDE       = "dde"
DISORDER  = "disorder"
DPC       = "dpc"
EAAC      = "eaac"
EGAAC     = "egaac"
GAAC      = "gaac"
GDPC      = "gdpc"
GTPC      = "gtpc"
GEARY     = "geary"
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
                        "gtpc", "gdpc", "dpc", "dde", "blosum62", "zscale"]

REST_ENCODINGS = ["aaindex", "psekraac"]

STRUC_ENCODINGS = ["disorder", "spinex", "psipred", "pssm"]

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
    """
    Gets the unique type, also referred to as ’subclass’, of an encoding. Ranges from 'easy'
    subclasses, such as lambda or window values to complex shapes, composed of several para-
    meters.
    :param encoding: The encoding.
    :param config: The snakemake config object.
    :return: The specific type of the given encoding.
    """
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
    """
    Gets the unique type of an encoding. The unique type is commonly the
    name of a parameter family, which is shared by all manifestations of
    a parametrized encoding.
    :param encoding: The encoding.
    :return: The unique type for a given encoding.
    """
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
    """
    Determines the type of input data to get the final datasets. Parameter-
    based encodings require more complex preceeding steps than encodings
    without parameters. Hence, the input depends on the encoding type.
    :param wildcards: The snakemake wildcards object.
    :param config: The snakemake config object.
    :return: The input data to compute the final (filtered) dataset.
    """
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


def get_encoding_description(encoding):
    descriptions = {
        AAC: "Amino Acid Composition",
        AAINDEX: "AAindex",
        APAAC: "Amphiphilic Pseudo-Amino Acid Composition",
        BINARY: "Binary",
        BLOSUM62: "Blosum62",
        CKSAAGP: "Composition of k-Spaced Amino Acid Group Pairs",
        CKSAAP: "Composition of k-spaced Amino Acid Pairs",
        CTDC: "Composition/Transition/Distribution - Composition",
        CTDD: "Composition/Transition/Distribution - Distribution",
        CTDT: "Composition/Transition/Distribution - Transition",
        CTRIAD: "Conjoint Triad",
        DDE: "Dipeptide Deviation from Expected Mean",
        DISORDER: "Disorder content and binary",
        DPC: "Di-Peptide Composition",
        EAAC: " Enhanced Amino Acid Composition",
        EGAAC: "Enhanced Grouped Amino Acid Composition",
        GAAC: "Grouped Amino Acid Composition",
        GDPC: "Grouped Di-Peptide Composition",
        GTPC: "Grouped Tri-Peptide Composition",
        GEARY: "Geary correlation",
        KSCTRIAD: "k-Spaced Conjoint Triad",
        MORAN: "Moran correlation",
        NMBROTO: "Normalized Moreau-Broto Autocorrelation",
        PAAC: "Pseudo-Amino Acid Composition",
        PSEKRAAC: "48 pseudo K-tuple reduced amino acids composition",
        PSIPRED: "PSIPRED - secondary structure content and binary",
        PSSM: "PSSM profile",
        QSORDER: "Quasi-sequence-order",
        SOCNUMBER: "Sequence-Order-Coupling Number",
        SPINEX: "SPINE-X - accessible solvent accessibility and torsion angle",
        TPC: "Tri-Peptide Composition",
        ZSCALE: "Z-Scale" }
    return descriptions[encoding]
