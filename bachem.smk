from utils.snakemake_config import WorkflowExecuter
import pandas as pd
from glob import glob

DATASET = "bachem"
# DATASETS = expand("bachem_window_length_{window_length}", window_length=range(6,21))
DATASETS = ["bachem_window_length_8"]
CORES = 4

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

rule all:
    input:
         expand("data/{normalized_dataset}/machine_learning/best_encodings.csv",
                normalized_dataset=DATASETS)

########################################################################################################################
############################################## DATASET CREATION ########################################################
########################################################################################################################

rule utils_sliding_windows:
    input:
         series_in=f"data/{DATASET}/series.yaml"
    output:
         fastas_out = "data/bachem_window_length_8/seqs.fasta",
         classes_out = "data/bachem_window_length_8/classes.txt"
         # fastas_out=expand("data/bachem_window_length_{window_length}/seqs.fasta", window_length=range(6,21)),
         # classes_out=expand("data/bachem_window_length_{window_length}/classes.txt", window_length=range(6,21))
    params:
         snakefile="nodes/utils/sliding_windows/Snakefile",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
############################################## PROFILE CREATION ########################################################
########################################################################################################################

rule util_multiple_sequence_alignment:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         classes_in="data/{normalized_dataset}/classes.txt"
    output:
         fasta_out="data/{normalized_dataset}/seqs_msa.fasta"
    params:
         snakefile="nodes/utils/multiple_sequence_alignment/Snakefile",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule util_secondary_structure_profile:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
         uniprot90_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                    "download_links/uniprot90_download_link.txt",
         psipred_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                  "download_links/psipred_download_link.txt",
         spineXpublic_download_link_in="nodes/utils/secondary_structure_profile/" + \
                                       "download_links/spineXpublic_download_link.txt",
         VSL2_download_link_in="nodes/utils/secondary_structure_profile/" + \
                               "download_links/VSL2_download_link.txt"
    output:
         fasta_anno_out="data/{normalized_dataset}/annotated_seqs.fasta",
         fasta_anno_msa_out="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         classes_anno="data/{normalized_dataset}/annotated_classes.txt",
         profiles_out=directory("data/{normalized_dataset}/profile/")
    params:
         snakefile="nodes/utils/secondary_structure_profile/Snakefile",
         configfile="nodes/utils/secondary_structure_profile/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule util_protein_structure_prediction:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
         download_link_in="nodes/utils/protein_structure_prediction/" + \
                          "download_links/raptorx_download_link.txt",
         license_key_in="nodes/utils/protein_structure_prediction/" + \
                        "download_links/modeller_license_key.txt"
    output:
         fasta_out="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_out="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdbs_out=directory("data/{normalized_dataset}/pdb/")
    params:
         snakefile="nodes/utils/protein_structure_prediction/Snakefile",
         configfile="nodes/utils/protein_structure_prediction/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
############################################## PARAMETER COMPUTATION ###################################################
########################################################################################################################

rule util_window_length:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta"
    output:
         length_out=expand("data/{{normalized_dataset}}/misc/{encoding}.yaml",
                           encoding=["ksctriad", "moran", "nmbroto", "geary",
                                     "qsorder", "socnumber", "eaac", "cksaagp",
                                     "cksaap", "apaac", "paac"])
    params:
         snakefile="nodes/utils/window_length/Snakefile",
         configfile="nodes/utils/window_length/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule util_dim_size:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta"
    output:
         length_out=expand("data/{{normalized_dataset}}/misc/ngram_{type}{size}.yaml",
                           type=["a","e","s"], size=[2,3])
    params:
         snakefile="nodes/utils/dim_size/Snakefile",
         configfile="nodes/utils/dim_size/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
############################################## MISC ENCODINGS ##########################################################
########################################################################################################################

rule encoding_pssm:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/pssm.csv"
    params:
         snakefile="nodes/encodings/pssm/Snakefile",
         configfile="nodes/encodings/pssm/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_zscale:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/zscale.csv"
    params:
         snakefile="nodes/encodings/zscale/Snakefile",
         configfile="nodes/encodings/zscale/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_tpc:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/tpc.csv"
    params:
         snakefile="nodes/encodings/tpc/Snakefile",
         configfile="nodes/encodings/tpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_gtpc:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/gtpc.csv"
    params:
         snakefile="nodes/encodings/gtpc/Snakefile",
         configfile="nodes/encodings/gtpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_gdpc:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/gdpc.csv"
    params:
         snakefile="nodes/encodings/gdpc/Snakefile",
         configfile="nodes/encodings/gdpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_gaac:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/gaac.csv"
    params:
         snakefile="nodes/encodings/gaac/Snakefile",
         configfile="nodes/encodings/gaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_egaac:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/egaac/egaac_window_{window_val}.csv",
                        window_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/egaac/Snakefile",
         configfile="nodes/encodings/egaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_dpc:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/dpc.csv"
    params:
         snakefile="nodes/encodings/dpc/Snakefile",
         configfile="nodes/encodings/dpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_dde:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/dde.csv"
    params:
         snakefile="nodes/encodings/dde/Snakefile",
         configfile="nodes/encodings/dde/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ctdt:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/ctdt.csv"
    params:
         snakefile="nodes/encodings/ctdt/Snakefile",
         configfile="nodes/encodings/ctdt/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ctdd:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/ctdd.csv"
    params:
         snakefile="nodes/encodings/ctdd/Snakefile",
         configfile="nodes/encodings/ctdd/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ctdc:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/ctdc.csv"
    params:
         snakefile="nodes/encodings/ctdc/Snakefile",
         configfile="nodes/encodings/ctdc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_blosum62:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/blosum62.csv"
    params:
         snakefile="nodes/encodings/blosum62/Snakefile",
         configfile="nodes/encodings/blosum62/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_binary:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/binary.csv"
    params:
         snakefile="nodes/encodings/binary/Snakefile",
         configfile="nodes/encodings/binary/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_aaindex:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/aaindex/aaindex_{aaindex}.csv",
                         aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/aaindex/Snakefile",
         configfile="nodes/encodings/aaindex/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_aac:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/aac.csv"
    params:
         snakefile="nodes/encodings/aac/Snakefile",
         configfile="nodes/encodings/aac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_fft:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         csv_in=expand("data/{{normalized_dataset}}/csv/aaindex/aaindex_{aaindex}.csv",
                       aaindex=get_aaindex())
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/fft/fft_{aaindex}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/fft/Snakefile",
         configfile="nodes/encodings/fft/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_waac:
    input:
         csv_in="data/{normalized_dataset}/csv/aac.csv",
         aaindex_in="apps/iFeature/data/AAindex.tsv"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/waac/waac_{aaindex}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/waac/Snakefile",
         configfile="nodes/encodings/waac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_flgc:
    input:
         csv_in="data/{normalized_dataset}/csv/aac.csv",
         aaindex_in="apps/iFeature/data/AAindex.tsv"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/flgc/flgc_{aaindex}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/flgc/Snakefile",
         configfile="nodes/encodings/flgc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_fldpc:
    input:
         csv_in="data/{normalized_dataset}/csv/dpc.csv",
         aaindex_in="apps/iFeature/data/AAindex.tsv"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/fldpc/fldpc_{aaindex}.csv",
                        aaindex=get_aaindex())
    params:
         snakefile="nodes/encodings/fldpc/Snakefile",
         configfile="nodes/encodings/fldpc/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ngram_a2:
    input:
         csv_in="data/{normalized_dataset}/csv/dpc.csv",
         length_in="data/{normalized_dataset}/misc/ngram_a2.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ngram_a2/ngram_a2_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         lsv_out=expand("data/{{normalized_dataset}}/csv/ngram_a2/ngram_a2_lsv_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         sv_out=expand("data/{{normalized_dataset}}/csv/ngram_a2/ngram_a2_sv_{dim}.csv",
                       dim=[1, 5, 20, 50, 100, 200, 300])
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ngram_a3:
    input:
         csv_in="data/{normalized_dataset}/csv/tpc.csv",
         length_in="data/{normalized_dataset}/misc/ngram_a3.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ngram_a3/ngram_a3_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         lsv_out=expand("data/{{normalized_dataset}}/csv/ngram_a3/ngram_a3_lsv_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         sv_out=expand("data/{{normalized_dataset}}/csv/ngram_a3/ngram_a3_sv_{dim}.csv",
                       dim=[1, 5, 20, 50, 100, 200, 300])
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ngram_e2:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/ngram_e2.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ngram_e2/ngram_e2_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         lsv_out=expand("data/{{normalized_dataset}}/csv/ngram_e2/ngram_e2_lsv_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         sv_out=expand("data/{{normalized_dataset}}/csv/ngram_e2/ngram_e2_sv_{dim}.csv",
                       dim=[1, 5, 20, 50, 100, 200, 300])
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ngram_e3:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/ngram_e3.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ngram_e3/ngram_e3_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         lsv_out=expand("data/{{normalized_dataset}}/csv/ngram_e3/ngram_e3_lsv_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         sv_out=expand("data/{{normalized_dataset}}/csv/ngram_e3/ngram_e3_sv_{dim}.csv",
                       dim=[1, 5, 20, 50, 100, 200, 300])
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ngram_s2:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/ngram_s2.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ngram_s2/ngram_s2_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         lsv_out=expand("data/{{normalized_dataset}}/csv/ngram_s2/ngram_s2_lsv_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         sv_out=expand("data/{{normalized_dataset}}/csv/ngram_s2/ngram_s2_sv_{dim}.csv",
                       dim=[1, 5, 20, 50, 100, 200, 300])
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ngram_s3:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/ngram_s3.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ngram_s3/ngram_s3_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         lsv_out=expand("data/{{normalized_dataset}}/csv/ngram_s3/ngram_s3_lsv_{dim}.csv",
                        dim=[1, 5, 20, 50, 100, 200, 300]),
         sv_out=expand("data/{{normalized_dataset}}/csv/ngram_s3/ngram_s3_sv_{dim}.csv",
                       dim=[1, 5, 20, 50, 100, 200, 300])
    params:
         snakefile="nodes/encodings/ngram/Snakefile",
         configfile="nodes/encodings/ngram/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_cgr:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/cgr/cgr_res_{{resolution}}_sf_{sfactor}.csv",
                        resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713])
    params:
         snakefile="nodes/encodings/cgr/Snakefile",
         configfile="nodes/encodings/cgr/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ctriad:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/ctriad.csv"
    params:
         snakefile="nodes/encodings/ctriad/Snakefile",
         configfile="nodes/encodings/ctriad/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_distance_frequency:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/distance_frequency/" + \
                        "distance_frequency_dn_{nterminal}_dc_{cterminal}.csv",
                        nterminal=[5, 10, 20, 50, 100],
                        cterminal=[5, 10, 20, 50, 100])
    params:
         snakefile="nodes/encodings/distance_frequency/Snakefile",
         configfile="nodes/encodings/distance_frequency/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_blomap:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out="data/{normalized_dataset}/csv/blomap.csv"
    params:
         snakefile="nodes/encodings/blomap/Snakefile",
         configfile="nodes/encodings/blomap/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
########################################### PARAM_BASED ENCODINGS ######################################################
########################################################################################################################

rule encoding_cksaagp:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/cksaagp.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/cksaagp/cksaagp_gap_{gap_val}.csv",
                        gap_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/cksaagp/Snakefile",
         configfile="nodes/encodings/cksaagp/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_socnumber:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/socnumber.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/socnumber/socnumber_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/socnumber/Snakefile",
         configfile="nodes/encodings/socnumber/config.yaml"
    resources:
         cores=4
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_qsorder:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/qsorder.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/qsorder/qsorder_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/qsorder/Snakefile",
         configfile="nodes/encodings/qsorder/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_nmbroto:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/nmbroto.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/nmbroto/Snakefile",
         configfile="nodes/encodings/nmbroto/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_moran:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/moran.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/moran/moran_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/moran/Snakefile",
         configfile="nodes/encodings/moran/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ksctriad:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/ksctriad.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/ksctriad/ksctriad_gap_{gap_val}.csv",
                        gap_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/ksctriad/Snakefile",
         configfile="nodes/encodings/ksctriad/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_geary:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/geary.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/geary/geary_nlag_{nlag_val}.csv",
                        nlag_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/geary/Snakefile",
         configfile="nodes/encodings/geary/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_eaac:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/eaac.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/eaac/eaac_window_{window_val}.csv",
                        window_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/eaac/Snakefile",
         configfile="nodes/encodings/eaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_cksaap:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/cksaap.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/cksaap/cksaap_gap_{gap_val}.csv",
                        gap_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/cksaap/Snakefile",
         configfile="nodes/encodings/cksaap/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_apaac:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/apaac.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/apaac/apaac_lambda_{lambda_val}.csv",
                        lambda_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/apaac/Snakefile",
         configfile="nodes/encodings/apaac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_paac:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         length_in="data/{normalized_dataset}/misc/paac.yaml"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/paac/paac_lambda_{lambda_val}.csv",
                        lambda_val=list(range(1, 31)))
    params:
         snakefile="nodes/encodings/paac/Snakefile",
         configfile="nodes/encodings/paac/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
############################################## PSEKRAAC ENCODINGS ######################################################
########################################################################################################################

rule encoding_psekraac_type16:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type16/psekraac_type16_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type16/Snakefile",
         configfile="nodes/encodings/psekraac_type16/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type15:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type15/psekraac_type15_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type15/Snakefile",
         configfile="nodes/encodings/psekraac_type15/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type14:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type14/psekraac_type14_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type14/Snakefile",
         configfile="nodes/encodings/psekraac_type14/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type13:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type13/psekraac_type13_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type13/Snakefile",
         configfile="nodes/encodings/psekraac_type13/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type12:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type12/psekraac_type12_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type12/Snakefile",
         configfile="nodes/encodings/psekraac_type12/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type11:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type11/psekraac_type11_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type11/Snakefile",
         configfile="nodes/encodings/psekraac_type11/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type10:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type10/psekraac_type10_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type10/Snakefile",
         configfile="nodes/encodings/psekraac_type10/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type9:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type9/psekraac_type9_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type9/Snakefile",
         configfile="nodes/encodings/psekraac_type9/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type8:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type8/psekraac_type8_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type8/Snakefile",
         configfile="nodes/encodings/psekraac_type8/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type7:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type7/psekraac_type7_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type7/Snakefile",
         configfile="nodes/encodings/psekraac_type7/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type6C:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type6C/psekraac_type6C_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type6C/Snakefile",
         configfile="nodes/encodings/psekraac_type6C/config.yaml"
    resources:
         cores=8
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type6B:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type6B/psekraac_type6B_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type6B/Snakefile",
         configfile="nodes/encodings/psekraac_type6B/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type6A:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type6A/psekraac_type6A_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type6A/Snakefile",
         configfile="nodes/encodings/psekraac_type6A/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type5:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type5/psekraac_type5_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type5/Snakefile",
         configfile="nodes/encodings/psekraac_type5/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type4:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type4/psekraac_type4_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type4/Snakefile",
         configfile="nodes/encodings/psekraac_type4/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type3B:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type3B/psekraac_type3B_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type3B/Snakefile",
         configfile="nodes/encodings/psekraac_type3B/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type3A:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type3A/psekraac_type3A_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type3A/Snakefile",
         configfile="nodes/encodings/psekraac_type3A/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type2:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type2/psekraac_type2_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type2/Snakefile",
         configfile="nodes/encodings/psekraac_type2/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_psekraac_type1:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/psekraac_type1/psekraac_type1_"
                        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/encodings/psekraac_type1/Snakefile",
         configfile="nodes/encodings/psekraac_type1/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
############################################ STRUCTURE-BASED ENCODINGS #################################################
########################################################################################################################

rule encoding_asa:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/asa.csv"
    params:
         snakefile="nodes/encodings/asa/Snakefile",
         configfile="nodes/encodings/asa/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ta:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/ta.csv"
    params:
         snakefile="nodes/encodings/ta/Snakefile",
         configfile="nodes/encodings/ta/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_ssec:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/ssec.csv"
    params:
         snakefile="nodes/encodings/ssec/Snakefile",
         configfile="nodes/encodings/ssec/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_sseb:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/sseb.csv"
    params:
         snakefile="nodes/encodings/sseb/Snakefile",
         configfile="nodes/encodings/sseb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_disorder:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/disorder.csv"
    params:
         snakefile="nodes/encodings/disorder/Snakefile",
         configfile="nodes/encodings/disorder/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_disorderb:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/disorderb.csv"
    params:
         snakefile="nodes/encodings/disorderb/Snakefile",
         configfile="nodes/encodings/disorderb/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_disorderc:
    input:
         fasta_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_classes.txt",
         profile="data/{normalized_dataset}/profile/"
    output:
         csv_out="data/{normalized_dataset}/csv/disorderc.csv"
    params:
         snakefile="nodes/encodings/disorderc/Snakefile",
         configfile="nodes/encodings/disorderc/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_qsar:
    input:
         fasta_in="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         csv_out="data/{normalized_dataset}/csv/qsar.csv"
    params:
         snakefile="nodes/encodings/qsar/Snakefile",
         configfile="nodes/encodings/qsar/config.yaml",
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_electrostatic_hull:
    input:
         fasta_in="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                        distance=[0,3,6,9,12])
    params:
         snakefile="nodes/encodings/electrostatic_hull/Snakefile",
         configfile="nodes/encodings/electrostatic_hull/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_distance_distribution:
    input:
         fasta_in="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         csv_out="data/{normalized_dataset}/csv/distance_distribution.csv"
    params:
         snakefile="nodes/encodings/distance_distribution/Snakefile",
         configfile="nodes/encodings/distance_distribution/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

rule encoding_delaunay:
    input:
         fasta_in="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_in="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         csv_out=expand("data/{{normalized_dataset}}/csv/delaunay/delaunay_{algorithm}.csv",
                        algorithm=["average_distance", "total_distance", "cartesian_product",
                                   "number_instances", "frequency_instances"])
    params:
         snakefile="nodes/encodings/delaunay/Snakefile",
         configfile="nodes/encodings/delaunay/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

########################################################################################################################
################################################ MACHINE LEARNING ######################################################
########################################################################################################################

rule machine_learning_top_encodings:
    input:
         # expand("data/{normalized_dataset}/csv/aaindex/aaindex_{aaindex}.csv",
         #        normalized_dataset=DATASETS,
         #        aaindex=get_aaindex()) +
         expand("data/{normalized_dataset}/csv/apaac/apaac_lambda_{lambda_val}.csv",
                normalized_dataset=DATASETS, lambda_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/cksaagp/cksaagp_gap_{gap_val}.csv",
         #        normalized_dataset=DATASETS, gap_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/socnumber/socnumber_nlag_{nlag_val}.csv",
         #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/qsorder/qsorder_nlag_{nlag_val}.csv",
         #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv",
         #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/moran/moran_nlag_{nlag_val}.csv",
         #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/ksctriad/ksctriad_gap_{gap_val}.csv",
         #        normalized_dataset=DATASETS, gap_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/geary/geary_nlag_{nlag_val}.csv",
         #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/eaac/eaac_window_{window_val}.csv",
         #        normalized_dataset=DATASETS, window_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/cksaap/cksaap_gap_{gap_val}.csv",
         #        normalized_dataset=DATASETS, gap_val=list(range(1, 31))) +
         # expand("data/{normalized_dataset}/csv/paac/paac_lambda_{lambda_val}.csv",
         #        normalized_dataset=DATASETS, lambda_val=list(range(1, 31))) +
         #
         # expand("data/{normalized_dataset}/csv/psekraac_type1/psekraac_type1_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type2/psekraac_type2_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type3A/psekraac_type3A_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type3B/psekraac_type3B_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type4/psekraac_type4_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type5/psekraac_type5_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type6A/psekraac_type6A_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type6B/psekraac_type6B_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type6C/psekraac_type6C_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type7/psekraac_type7_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type8/psekraac_type8_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type9/psekraac_type9_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type10/psekraac_type10_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type11/psekraac_type11_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type12/psekraac_type12_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type13/psekraac_type13_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type14/psekraac_type14_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type15/psekraac_type15_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         # expand("data/{normalized_dataset}/csv/psekraac_type16/psekraac_type16_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        normalized_dataset=DATASETS,
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))) +
         #
         # expand("data/{normalized_dataset}/csv/fft/fft_{aaindex}.csv",
         #        normalized_dataset=DATASETS, aaindex=get_aaindex()) +
         #
         # expand("data/{normalized_dataset}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
         #        normalized_dataset=DATASETS,
         #        resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]) +
         #
         # expand("data/{normalized_dataset}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
         #        normalized_dataset=DATASETS,
         #        distance=[0,3,6,9,12]) +
         # expand("data/{normalized_dataset}/csv/delaunay/delaunay_{algorithm}.csv",
         #        normalized_dataset=DATASETS,
         #        algorithm=["average_distance", "total_distance", "cartesian_product",
         #                   "number_instances", "frequency_instances"]) +
         # expand("data/{normalized_dataset}/csv/waac/waac_{aaindex}.csv",
         #        normalized_dataset=DATASETS,
         #        aaindex=get_aaindex()) +
         # expand("data/{normalized_dataset}/csv/flgc/flgc_{aaindex}.csv",
         #        normalized_dataset=DATASETS,
         #        aaindex=get_aaindex()) +
         # expand("data/{normalized_dataset}/csv/fldpc/fldpc_{aaindex}.csv",
         #        normalized_dataset=DATASETS,
         #        aaindex=get_aaindex()) +
         #
         # expand("data/{normalized_dataset}/csv/ngram_a2/ngram_a2_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_a2/ngram_a2_lsv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_a2/ngram_a2_sv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_a3/ngram_a3_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_a3/ngram_a3_lsv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_a3/ngram_a3_sv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_e2/ngram_e2_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_e2/ngram_e2_lsv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_e2/ngram_e2_sv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_e3/ngram_e3_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_e3/ngram_e3_lsv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_e3/ngram_e3_sv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_s2/ngram_s2_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_s2/ngram_s2_lsv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_s2/ngram_s2_sv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_s3/ngram_s3_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_s3/ngram_s3_lsv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         # expand("data/{normalized_dataset}/csv/ngram_s3/ngram_s3_sv_{dim}.csv",
         #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300]) +
         #
         # expand("data/{normalized_dataset}/csv/disorder.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/disorderb.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/disorderc.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/aac.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/binary.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/blosum62.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/ctdc.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/ctdd.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/ctdt.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/ctriad.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/dde.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/dpc.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/gaac.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/gdpc.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/gtpc.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/tpc.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/zscale.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/pssm.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/sseb.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/ssec.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/ta.csv",
         #        normalized_dataset=DATASETS) +
         # expand("data/{normalized_dataset}/csv/asa.csv",
         #        normalized_dataset=DATASETS) +
         expand("data/{normalized_dataset}/csv/blomap.csv",
                normalized_dataset=DATASETS) +
         expand("data/{normalized_dataset}/csv/distance_distribution.csv",
                normalized_dataset=DATASETS)
         # expand("data/{normalized_dataset}/csv/qsar.csv",
         #       normalized_dataset=DATASETS)
    output:
         csv_out="data/{normalized_dataset}/machine_learning/top_encodings.csv"
    params:
         snakefile="nodes/encodings/distance_distribution/Snakefile",
         configfile="nodes/encodings/distance_distribution/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile):
             shell(f"""snakemake -s {{params.snakefile}} --cores {CORES} --configfile {{params.configfile}}""")

# https://scikit-learn.org/stable/modules/ensemble.html#voting-classifier