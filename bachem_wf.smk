from modlamp.core import read_fasta
from utils.snakemake_config import WorkflowExecuter
import pandas as pd
import yaml
import sys

DATASET = config["dataset"]
CORES = 8

include:
    "maximum_window_length.smk"

include: 
    "param_based_encodings.smk"

include:
    "misc_encodings.smk"

include:
    "structure_based_encodings.smk"

include:
    "psekraac_encodings.smk"

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

max_vals = {}
for encoding in ["ksctriad", "moran", "nmbroto", "geary","qsorder",
                 "socnumber", "eaac", "cksaagp", "cksaap", "apaac", "paac"]:
    try:
        with open(f"data/{DATASET}/misc/{encoding}.yaml") as f:
            max_vals[encoding] = yaml.safe_load(f) + 1  # range is exclusive
    except FileNotFoundError:
        sys.exit("""
        Please run node window_length beforehand (or set them manually).
        See, e.g., apps/iFeature/codes/KSCTriad.py for details.
        """)

rule all:
    input:
        expand(f"data/{DATASET}/csv/aaindex/aaindex_{{aaindex}}.csv", aaindex=get_aaindex()),
        expand(f"data/{DATASET}/csv/apaac/apaac_lambda_{{lambda_val}}.csv", lambda_val=list(range(1, max_vals["apaac"]))),
        expand(f"data/{DATASET}/csv/cksaagp/cksaagp_gap_{{gap_val}}.csv", gap_val=list(range(1, max_vals["cksaagp"]))),
        expand(f"data/{DATASET}/csv/cksaap/cksaap_gap_{{gap_val}}.csv", gap_val=list(range(1, max_vals["cksaap"]))),
        expand(f"data/{DATASET}/csv/eaac/eaac_window_{{window_val}}.csv", window_val=list(range(1, max_vals["eaac"]))),
        expand(f"data/{DATASET}/csv/egaac/egaac_window_{{window_val}}.csv", window_val=list(range(1, 31))),
        expand(f"data/{DATASET}/csv/geary/geary_nlag_{{nlag_val}}.csv", nlag_val=list(range(1, max_vals["geary"]))),
        expand(f"data/{DATASET}/csv/ksctriad/ksctriad_gap_{{gap_val}}.csv", gap_val=list(range(1, max_vals["ksctriad"]))),
        expand(f"data/{DATASET}/csv/moran/moran_nlag_{{nlag_val}}.csv", nlag_val=list(range(1, max_vals["moran"]))),
        expand(f"data/{DATASET}/csv/nmbroto/nmbroto_nlag_{{nlag_val}}.csv", nlag_val=list(range(1, max_vals["nmbroto"]))),
        expand(f"data/{DATASET}/csv/paac/paac_lambda_{{lambda_val}}.csv", lambda_val=list(range(1, max_vals["paac"]))),
        expand(f"data/{DATASET}/csv/qsorder/qsorder_nlag_{{nlag_val}}.csv", nlag_val=list(range(1, max_vals["qsorder"]))),
        expand(f"data/{DATASET}/csv/socnumber/socnumber_nlag_{{nlag_val}}.csv", nlag_val=list(range(1, max_vals["socnumber"]))),
        f"data/{DATASET}/csv/disorder.csv",
        f"data/{DATASET}/csv/disorderb.csv",
        f"data/{DATASET}/csv/disorderc.csv",
        f"data/{DATASET}/csv/aac.csv",
        f"data/{DATASET}/csv/binary.csv",
        f"data/{DATASET}/csv/blosum62.csv",
        f"data/{DATASET}/csv/ctdc.csv",
        f"data/{DATASET}/csv/ctdd.csv",
        f"data/{DATASET}/csv/ctdt.csv",
        f"data/{DATASET}/csv/ctriad.csv",
        f"data/{DATASET}/csv/dde.csv",
        f"data/{DATASET}/csv/dpc.csv",
        f"data/{DATASET}/csv/gaac.csv",
        f"data/{DATASET}/csv/gdpc.csv",
        f"data/{DATASET}/csv/gtpc.csv",
        f"data/{DATASET}/csv/tpc.csv",
        f"data/{DATASET}/csv/zscale.csv",
        f"data/{DATASET}/csv/pssm.csv",
        f"data/{DATASET}/csv/sseb.csv",
        f"data/{DATASET}/csv/ssec.csv",
        f"data/{DATASET}/csv/ta.csv",
        f"data/{DATASET}/csv/asa.csv",
        expand(f"data/{DATASET}/csv/psekraac_type1/psekraac_type1_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type2/psekraac_type2_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type3A/psekraac_type3A_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type3B/psekraac_type3B_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type4/psekraac_type4_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type5/psekraac_type5_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type6A/psekraac_type6A_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type6B/psekraac_type6B_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type6C/psekraac_type6C_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type7/psekraac_type7_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type8/psekraac_type8_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type9/psekraac_type9_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type10/psekraac_type10_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type11/psekraac_type11_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type12/psekraac_type12_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type13/psekraac_type13_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type14/psekraac_type14_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type15/psekraac_type15_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/psekraac_type16/psekraac_type16_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
        expand(f"data/{DATASET}/csv/fft/fft_{{aaindex}}.csv", aaindex=get_aaindex()),
        expand(f"data/{DATASET}/csv/cgr/cgr_res_{{resolution}}_sf_{{sfactor}}.csv",
               resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),
        f"data/{DATASET}/csv/distance_distribution.csv",
        f"data/{DATASET}/csv/qsar.csv",
        expand(f"data/{DATASET}/csv/electrostatic_hull/electrostatic_hull_{{distance}}.csv",
               distance=[0,3,6,9,12])

rule util_multiple_sequence_alignment:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt"
    output:
         fasta_out=f"data/{DATASET}/seqs_msa.fasta"
    params:
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} {{output.fasta_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule util_secondary_structure_profile:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         fasta_msa_in=f"data/{DATASET}/seqs_msa.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         uniprot90_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                    "download_links/uniprot90_download_link.txt",
         psipred_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                  "download_links/psipred_download_link.txt",
         spineXpublic_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                       "download_links/spineXpublic_download_link.txt",
         VSL2_download_link_in="nodes/utils/secondary_structure_profile/" + \
                               "download_links/VSL2_download_link.txt"
    output:
         fasta_anno_out=f"data/{DATASET}/annotated_seqs.fasta",
         fasta_anno_msa_out=f"data/{DATASET}/annotated_seqs_msa.fasta",
         classes_anno=f"data/{DATASET}/annotated_classes.txt",
         profiles_out=expand(f"data/{DATASET}/profile/{{seq_name}}.{{ftype}}",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
                         ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"])
    params:
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} \
                            {{output.fasta_anno_out}} {{output.fasta_anno_msa_out}} {{output.classes_anno}} {{output.profiles_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")

rule util_protein_structure_prediction:
    input:
         fasta_in=f"data/{DATASET}/seqs.fasta",
         classes_in=f"data/{DATASET}/classes.txt",
         download_link_in="nodes/utils/protein_structure_prediction/" + \
                          "download_links/raptorx_download_link.txt",
         license_key_in="nodes/utils/protein_structure_prediction/" + \
                        "download_links/modeller_license_key.txt"
    output:
         fasta_out=f"data/{DATASET}/annotated_pdbs_seqs.fasta",
         classes_out=f"data/{DATASET}/annotated_pdbs_classes.txt",
         pdbs_out=expand(f"data/{DATASET}/pdb/{{seq_name}}.pdb",
                         seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1])
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} \
                            {{output.fasta_out}} {{output.classes_out}} {{output.pdbs_out}} \
                            --cores {CORES} \
                            --directory $PWD \
                            --configfile {{params.configfile}}""")



