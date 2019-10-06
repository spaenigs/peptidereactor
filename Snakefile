import os
import pandas as pd

config["dataset"] = "neuropeptides_ds3"
config["global_workdir"] = os.getcwd() + "/"

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()[:5]

rule all:
    input:
        "data/neuropeptides_ds3/csv/disorder.csv",
        "data/neuropeptides_ds3/csv/disorderb.csv",
        "data/neuropeptides_ds3/csv/disorderc.csv",
        "data/neuropeptides_ds3/csv/aac.csv",
        expand("data/neuropeptides_ds3/csv/aaindex/aaindex_{aaindex}.csv", aaindex=get_aaindex()),
        expand("data/neuropeptides_ds3/csv/apaac/apaac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 4))),
        "data/neuropeptides_ds3/csv/binary.csv",
        "data/neuropeptides_ds3/csv/blosum62.csv",
        expand("data/neuropeptides_ds3/csv/cksaagp/cksaagp_gap_{gap_val}.csv", gap_val=list(range(1, 4))),
        expand("data/neuropeptides_ds3/csv/cksaap/cksaap_gap_{gap_val}.csv", gap_val=list(range(1, 4))),
        "data/neuropeptides_ds3/csv/ctdc.csv",
        "data/neuropeptides_ds3/csv/ctdd.csv",
        "data/neuropeptides_ds3/csv/ctdt.csv",
        "data/neuropeptides_ds3/csv/ctriad.csv",
        "data/neuropeptides_ds3/csv/dde.csv",
        "data/neuropeptides_ds3/csv/dpc.csv",
        expand("data/neuropeptides_ds3/csv/eaac/eaac_window_{window_val}.csv", window_val=list(range(1, 4))),
        expand("data/neuropeptides_ds3/csv/egaac/egaac_window_{window_val}.csv", window_val=list(range(1, 4))),
        "data/neuropeptides_ds3/csv/gaac.csv",
        "data/neuropeptides_ds3/csv/gdpc.csv",
        expand("data/neuropeptides_ds3/csv/geary/geary_nlag_{nlag_val}.csv", nlag_val=list(range(1, 4))),
        "data/neuropeptides_ds3/csv/gtpc.csv",
        expand("data/neuropeptides_ds3/csv/ksctriad/ksctriad_gap_{gap_val}.csv", gap_val=list(range(1, 4))),
        expand("data/neuropeptides_ds3/csv/moran/moran_nlag_{nlag_val}.csv", nlag_val=list(range(1, 4))),
        expand("data/neuropeptides_ds3/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv", nlag_val=list(range(1, 4))),
        expand("data/neuropeptides_ds3/csv/paac/paac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 4))),
        expand("data/neuropeptides_ds3/csv/qsorder/qsorder_nlag_{nlag_val}.csv", nlag_val=list(range(1, 4))),
        "data/neuropeptides_ds3/csv/tpc.csv",
        "data/neuropeptides_ds3/csv/zscale.csv",
        expand("data/neuropeptides_ds3/csv/psekraac/type1/psekraac_type1_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["lambda-correlation"], raac_val=[2],
               ktuple_val=[2], lambda_val=[2]),
               # sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               # ktuple_val=list(range(1,4)), lambda_val=list(range(1,7)))
        "data/neuropeptides_ds3/csv/pssm.csv",
        expand("data/neuropeptides_ds3/csv/socnumber/socnumber_nlag_{nlag_val}.csv", nlag_val=list(range(1, 4))),
        "data/neuropeptides_ds3/csv/sseb.csv"

rule encoding_sseb:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/sseb.csv"
    params:
         subworkflow="sseb",
         snakefile="nodes/encodings/sseb/Snakefile",
         configfile="nodes/encodings/sseb/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_socnumber:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/socnumber/socnumber_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 4)))
    params:
         subworkflow="socnumber",
         snakefile="nodes/encodings/socnumber/Snakefile",
         configfile="nodes/encodings/socnumber/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_pssm:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/pssm.csv"
    params:
         subworkflow="pssm",
         snakefile="nodes/encodings/pssm/Snakefile",
         configfile="nodes/encodings/pssm/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_psekraac_type_1:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/psekraac/type1/psekraac_type1_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["lambda-correlation"], raac_val=[2],
                        ktuple_val=[2], lambda_val=[2])
                        # sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        # ktuple_val=list(range(1,4)), lambda_val=list(range(1,7)))
    params:
         subworkflow="psekraac_type1",
         snakefile="nodes/encodings/psekraac_type1/Snakefile",
         configfile="nodes/encodings/psekraac_type1/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"


rule encoding_zscale:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/zscale.csv"
    params:
         subworkflow="zscale",
         snakefile="nodes/encodings/zscale/Snakefile",
         configfile="nodes/encodings/zscale/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_tpc:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/tpc.csv"
    params:
         subworkflow="tpc",
         snakefile="nodes/encodings/tpc/Snakefile",
         configfile="nodes/encodings/tpc/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_qsorder:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/qsorder/qsorder_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 4)))
    params:
         subworkflow="qsorder",
         snakefile="nodes/encodings/qsorder/Snakefile",
         configfile="nodes/encodings/qsorder/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_paac:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/paac/paac_lambda_{lambda_val}.csv",
                        lambda_val=list(range(1, 4)))
    params:
         subworkflow="paac",
         snakefile="nodes/encodings/paac/Snakefile",
         configfile="nodes/encodings/paac/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_nmbroto:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 4)))
    params:
         subworkflow="moran",
         snakefile="nodes/encodings/nmbroto/Snakefile",
         configfile="nodes/encodings/nmbroto/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_moran:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/moran/moran_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 4)))
    params:
         subworkflow="moran",
         snakefile="nodes/encodings/moran/Snakefile",
         configfile="nodes/encodings/moran/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_ksctriad:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/ksctriad/ksctriad_gap_{gap_val}.csv",
                        gap_val=list(range(1, 4)))
    params:
         subworkflow="ksctriad",
         snakefile="nodes/encodings/ksctriad/Snakefile",
         configfile="nodes/encodings/ksctriad/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_gtpc:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/gtpc.csv"
    params:
         subworkflow="gtpc",
         snakefile="nodes/encodings/gtpc/Snakefile",
         configfile="nodes/encodings/gtpc/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_geary:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/geary/geary_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 4)))
    params:
         subworkflow="geary",
         snakefile="nodes/encodings/geary/Snakefile",
         configfile="nodes/encodings/geary/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_gdpc:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/gdpc.csv"
    params:
         subworkflow="gdpc",
         snakefile="nodes/encodings/gdpc/Snakefile",
         configfile="nodes/encodings/gdpc/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_gaac:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/gaac.csv"
    params:
         subworkflow="gaac",
         snakefile="nodes/encodings/gaac/Snakefile",
         configfile="nodes/encodings/gaac/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_egaac:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/egaac/egaac_window_{window_val}.csv",
                        window_val=list(range(1, 4)))
    params:
         subworkflow="egaac",
         snakefile="nodes/encodings/egaac/Snakefile",
         configfile="nodes/encodings/egaac/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_eaac:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/eaac/eaac_window_{window_val}.csv",
                        window_val=list(range(1, 4)))
    params:
         subworkflow="eaac",
         snakefile="nodes/encodings/eaac/Snakefile",
         configfile="nodes/encodings/eaac/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_dpc:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/dpc.csv"
    params:
         subworkflow="dpc",
         snakefile="nodes/encodings/dpc/Snakefile",
         configfile="nodes/encodings/dpc/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_dde:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/dde.csv"
    params:
         subworkflow="dde",
         snakefile="nodes/encodings/dde/Snakefile",
         configfile="nodes/encodings/dde/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_ctriad:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/ctriad.csv"
    params:
         subworkflow="ctriad",
         snakefile="nodes/encodings/ctriad/Snakefile",
         configfile="nodes/encodings/ctriad/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_ctdt:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/ctdt.csv"
    params:
         subworkflow="ctdt",
         snakefile="nodes/encodings/ctdt/Snakefile",
         configfile="nodes/encodings/ctdt/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_ctdd:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/ctdd.csv"
    params:
         subworkflow="ctdd",
         snakefile="nodes/encodings/ctdd/Snakefile",
         configfile="nodes/encodings/ctdd/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_ctdc:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/ctdc.csv"
    params:
         subworkflow="ctdc",
         snakefile="nodes/encodings/ctdc/Snakefile",
         configfile="nodes/encodings/ctdc/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_cksaap:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/cksaap/cksaap_gap_{gap_val}.csv",
                        gap_val=list(range(1, 4)))
    params:
         subworkflow="cksaap",
         snakefile="nodes/encodings/cksaap/Snakefile",
         configfile="nodes/encodings/cksaap/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_cksaagp:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/cksaagp/cksaagp_gap_{gap_val}.csv",
                        gap_val=list(range(1, 4)))
    params:
         subworkflow="cksaagp",
         snakefile="nodes/encodings/cksaagp/Snakefile",
         configfile="nodes/encodings/cksaagp/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_blosum62:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/blosum62.csv"
    params:
         subworkflow="blosum62",
         snakefile="nodes/encodings/blosum62/Snakefile",
         configfile="nodes/encodings/blosum62/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_binary:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs_msa.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/binary.csv"
    params:
         subworkflow="binary",
         snakefile="nodes/encodings/binary/Snakefile",
         configfile="nodes/encodings/binary/config.yaml"
    script:
         "utils/subworkflow.py"

rule encoding_apaac:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/apaac/apaac_lambda_{lambda_val}.csv",
                        lambda_val=list(range(1, 4)))
    params:
         subworkflow="apaac",
         snakefile="nodes/encodings/apaac/Snakefile",
         configfile="nodes/encodings/apaac/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_aaindex:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out=expand("data/neuropeptides_ds3/csv/aaindex/aaindex_{aaindex}.csv",
                         aaindex=get_aaindex())
    params:
         subworkflow="aaindex",
         snakefile="nodes/encodings/aaindex/Snakefile",
         configfile="nodes/encodings/aaindex/config.yaml"
    resources:
         cores=4
    script:
         "utils/subworkflow.py"

rule encoding_disorder:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorder.csv"
    params:
         subworkflow="disorder",
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_disorderb:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs_msa.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorderb.csv"
    params:
         subworkflow="disorderb",
         snakefile="nodes/encodings/disorderb/Snakefile",
         configfile="nodes/encodings/disorderb/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_disorderc:
    input:
         fasta_in="data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt",
         profile=f"data/{config['dataset']}/profile"
    output:
         csv_out="data/neuropeptides_ds3/csv/disorderc.csv"
    params:
         subworkflow="disorderc",
         snakefile="nodes/encodings/disorderc/Snakefile",
         configfile="nodes/encodings/disorderc/config.yaml",
    script:
         "utils/subworkflow.py"

rule encoding_aac:
    input:
         fasta_in= "data/neuropeptides_ds3/annotated_seqs.fasta",
         classes_in="data/neuropeptides_ds3/annotated_classes.txt"
    output:
         csv_out="data/neuropeptides_ds3/csv/aac.csv"
    params:
         subworkflow="aac",
         snakefile="nodes/encodings/aac/Snakefile",
         configfile="nodes/encodings/aac/config.yaml"
    script:
         "utils/subworkflow.py"
