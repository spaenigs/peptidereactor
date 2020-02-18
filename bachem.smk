from proteinreactor.workflow_executer import WorkflowExecuter, MetaWorkflowExecuter
import pandas as pd

WINDOW_LENGTHS = [8,11,13,15,17,20]
DATASETS = \
    ["protease_window_length_8", "protease_window_length_8_complete"] + \
    expand(["bachem_window_length_{window_length}", "bachem_window_length_{window_length}_complete"],
           window_length=WINDOW_LENGTHS)

# WINDOW_LENGTHS = [13]#[8,11,15,20]
# DATASETS = ["protease_window_length_8", "bachem_window_length_13"]

CORES = int(config["cores"])

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

rule all:
    input:
         f"data/bachem/plots/filtered_datasets.png",
         f"data/protease/plots/filtered_datasets.png",
         # "data/protease/machine_learning/top_encodings.csv",
         # "data/bachem/machine_learning/top_encodings.csv",
         # "data/protease/machine_learning/ensemble_results_validated.csv",
         # "data/protease/machine_learning/ensemble_results_validated_tuned_hp.csv"

########################################################################################################################
############################################## DATASET CREATION ########################################################
########################################################################################################################

rule utils_sliding_windows:
    input:
         series_in=f"data/bachem/series.yaml"
    output:
         fastas_out=expand(f"data/bachem_window_length_{{window_length}}/seqs.fasta", window_length=WINDOW_LENGTHS),
         classes_out=expand(f"data/bachem_window_length_{{window_length}}/classes.txt", window_length=WINDOW_LENGTHS),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="bachem") as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_sliding_windows_complete:
    input:
         series_in=f"data/bachem/series.yaml"
    output:
         fastas_out=expand(f"data/bachem_window_length_{{window_length}}_complete/seqs.fasta", window_length=WINDOW_LENGTHS),
         classes_out=expand(f"data/bachem_window_length_{{window_length}}_complete/classes.yaml", window_length=WINDOW_LENGTHS),
         classes_idx_out=expand(f"data/bachem_window_length_{{window_length}}_complete/classes.txt", window_length=WINDOW_LENGTHS),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows_complete.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="bachem") as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_protein_dataset_creation:
    input:
         dataset_in="data/protease/impensData.txt",
    output:
         fasta_out="data/protease_window_length_8/seqs.fasta",
         classes_out="data/protease_window_length_8/classes.txt"
    params:
         snakefile="nodes/utils/protein_dataset_creation/protein_dataset_creation.smk",
         configfile="nodes/utils/protein_dataset_creation/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_protein_dataset_creation_complete:
    input:
         dataset_in="data/protease/impensData.txt",
         ids_file_in="data/protease/impens_ids.txt"
    output:
         fasta_out="data/protease_window_length_8_complete/seqs_orig.fasta",
         classes_out="data/protease_window_length_8_complete/classes_orig.txt",
         fasta_complete_out="data/protease_window_length_8_complete/seqs.fasta",
         classes_yaml_out="data/protease_window_length_8_complete/classes.yaml",
         classes_idx_out="data/protease_window_length_8_complete/classes.txt"
    params:
         snakefile="nodes/utils/protein_dataset_creation/protein_dataset_creation_complete.smk",
         configfile="nodes/utils/protein_dataset_creation/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

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
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="bachem") as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

########################################################################################################################
############################################ SEQUENCE-BASED ENCODINGS #################################################
########################################################################################################################

rule meta_workflow_sequence_based_encodings:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
    output:
         zscale_out="data/{normalized_dataset,.*?\d+}/csv/zscale.csv",
         tpc_out="data/{normalized_dataset,.*?\d+}/csv/tpc.csv",
         gtpc_out="data/{normalized_dataset,.*?\d+}/csv/gtpc.csv",
         gdpc_out="data/{normalized_dataset,.*?\d+}/csv/gdpc.csv",
         gaac_out="data/{normalized_dataset,.*?\d+}/csv/gaac.csv",
         dpc_out="data/{normalized_dataset,.*?\d+}/csv/dpc.csv",
         dde_out="data/{normalized_dataset,.*?\d+}/csv/dde.csv",
         ctdt_out="data/{normalized_dataset,.*?\d+}/csv/ctdt.csv",
         ctdd_out="data/{normalized_dataset,.*?\d+}/csv/ctdd.csv",
         ctdc_out="data/{normalized_dataset,.*?\d+}/csv/ctdc.csv",
         blosum62_out="data/{normalized_dataset,.*?\d+}/csv/blosum62.csv",
         binary_out="data/{normalized_dataset,.*?\d+}/csv/binary.csv",
         aac_out="data/{normalized_dataset,.*?\d+}/csv/aac.csv",
         ctriad_out="data/{normalized_dataset,.*?\d+}/csv/ctriad.csv",
         blomap_out="data/{normalized_dataset,.*?\d+}/csv/blomap.csv",
         egaac_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/egaac/egaac_window_{window_val}.csv", window_val=list(range(1, 31))),
         aaindex_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/aaindex/aaindex_{aaindex}.csv", aaindex=get_aaindex()),
         fft_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/fft/fft_{aaindex}.csv", aaindex=get_aaindex()),
         waac_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/waac/waac_{aaindex}.csv", aaindex=get_aaindex()),
         flgc_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/flgc/flgc_{aaindex}.csv", aaindex=get_aaindex()),
         fldpc_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/fldpc/fldpc_{aaindex}.csv", aaindex=get_aaindex()),
         ngram_a2_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_a2/ngram_a2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_a2_lsv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_a2/ngram_a2_lsv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_a2_sv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_a2/ngram_a2_sv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_a3_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_a3/ngram_a3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_a3_lsv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_a3/ngram_a3_lsv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_a3_sv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_a3/ngram_a3_sv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e2_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_e2/ngram_e2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e2_lsv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_e2/ngram_e2_lsv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e2_sv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_e2/ngram_e2_sv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e3_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_e3/ngram_e3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e3_lsv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_e3/ngram_e3_lsv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e3_sv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_e3/ngram_e3_sv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s2_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_s2/ngram_s2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s2_lsv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_s2/ngram_s2_lsv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s2_sv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_s2/ngram_s2_sv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s3_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_s3/ngram_s3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s3_lsv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_s3/ngram_s3_lsv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s3_sv_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ngram_s3/ngram_s3_sv_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         cgr_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
                     resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),
         distance_frequency_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/distance_frequency/distance_frequency_dn_{nterminal}_dc_{cterminal}.csv",
                     nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),
         cksaap_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/cksaap/cksaap_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         cksaagp_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/cksaagp/cksaagp_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         socnumber_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/socnumber/socnumber_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         qsorder_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/qsorder/qsorder_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         nmbroto_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         moran_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/moran/moran_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         ksctriad_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/ksctriad/ksctriad_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         geary_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/geary/geary_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         eaac_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/eaac/eaac_window_{window_val}.csv", window_val=list(range(1, 31))),
         cksaap=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/cksaap/cksaap_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         apaac_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/apaac/apaac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 31))),
         paac_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/paac/paac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 31))),
         psekraac_type16_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type16/psekraac_type16_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type15_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type15/psekraac_type15_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type14_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type14/psekraac_type14_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type13_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type13/psekraac_type13_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type12_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type12/psekraac_type12_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type11_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type11/psekraac_type11_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type10_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type10/psekraac_type10_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type9_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type9/psekraac_type9_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type8_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type8/psekraac_type8_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type7_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type7/psekraac_type7_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type6C_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type6C/psekraac_type6C_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type6B_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type6B/psekraac_type6B_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type6A_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type6A/psekraac_type6A_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type5_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type5/psekraac_type5_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type4_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type4/psekraac_type4_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type3B_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type3B/psekraac_type3B_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type3A_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type3A/psekraac_type3A_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type2_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type2/psekraac_type2_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type1_out=\
              expand("data/{{normalized_dataset,.*?\d+}}/csv/psekraac_type1/psekraac_type1_"
                     "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
                     sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
                     ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/meta_workflows/sequence_based_encodings/Snakefile",
         configfile="nodes/meta_workflows/sequence_based_encodings/config.yaml"
    priority:
         900
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

########################################################################################################################
############################################ STRUCTURE-BASED ENCODINGS #################################################
########################################################################################################################

rule meta_workflow_structure_based_profile:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt"
    output:
         fasta_anno_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_seqs.fasta",
         classes_anno="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?(1[2-9]|2\d)}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?(1[2-9]|2\d)}/pdb/")
    params:
         snakefile="nodes/meta_workflows/structure_based_profile/structure_based_profile.smk",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    priority:
         1000
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule meta_workflow_structure_based_encodings:
    input:
         fasta_anno="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_anno="data/{normalized_dataset}/annotated_classes.txt",
         fasta_anno_msa="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         profile_dir="data/{normalized_dataset}/profile/",
         fasta_anno_pdbs="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         pssm_out="data/{normalized_dataset}/csv/pssm.csv",
         asa_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/asa.csv",
         ta_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ta.csv",
         ssec_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ssec.csv",
         sseb_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/sseb.csv",
         disorder_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/disorder.csv",
         disorderb_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/disorderb.csv",
         disorderc_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/disorderc.csv",
         qsar_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/qsar.csv",
         electrostatic_hull_out=\
              expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                     distance=[0,3,6,9,12]),
         distance_distribution_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/distance_distribution.csv",
         delaunay_out=\
              expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/delaunay/delaunay_{algorithm}.csv",
                     algorithm=["average_distance", "total_distance", "cartesian_product",
                                "number_instances", "frequency_instances"])
    params:
         snakefile="nodes/meta_workflows/structure_based_encodings/structure_based_encodings.smk",
         configfile="nodes/meta_workflows/structure_based_encodings/config.yaml"
    priority:
         800
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule meta_workflow_structure_based_profile_windowed:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_idx_in="data/{normalized_dataset}/classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml"
    output:
         fasta_anno_out="data/{normalized_dataset,.*?[a-z]}/annotated_seqs.fasta",
         classes_anno_idx_out="data/{normalized_dataset,.*?[a-z]}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?[a-z]}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?[a-z]}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?[a-z]}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_idx_out="data/{normalized_dataset,.*?[a-z]}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?[a-z]}/pdb/")
    params:
         snakefile="nodes/meta_workflows/structure_based_profile/structure_based_profile_windowed.smk",
         configfile="nodes/utils/multiple_sequence_alignment/config.yaml"
    priority:
         1000
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule meta_workflow_structure_based_encodings_windowed:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_idx_in="data/{normalized_dataset}/classes.txt",
         classes_in="data/{normalized_dataset}/classes.yaml",
         fasta_anno="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_anno_idx="data/{normalized_dataset}/annotated_classes.txt",
         fasta_anno_msa="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         profile_dir="data/{normalized_dataset}/profile/",
         fasta_anno_pdbs="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_idx="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         pdb_dir="data/{normalized_dataset}/pdb/"
    output:
         asa_out="data/{normalized_dataset,.*?[a-z]}/csv/asa.csv",
         ta_out="data/{normalized_dataset,.*?[a-z]}/csv/ta.csv",
         ssec_out="data/{normalized_dataset,.*?[a-z]}/csv/ssec.csv",
         sseb_out="data/{normalized_dataset,.*?[a-z]}/csv/sseb.csv",
         disorder_out="data/{normalized_dataset,.*?[a-z]}/csv/disorder.csv",
         disorderb_out="data/{normalized_dataset,.*?[a-z]}/csv/disorderb.csv",
         disorderc_out="data/{normalized_dataset,.*?[a-z]}/csv/disorderc.csv",
         qsar_out="data/{normalized_dataset,.*?[a-z]}/csv/qsar.csv",
         electrostatic_hull_out=\
              expand("data/{{normalized_dataset,.*?[a-z]}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                     distance=[0,3,6,9,12]),
         distance_distribution_out="data/{normalized_dataset,.*?[a-z]}/csv/distance_distribution.csv",
         delaunay_out=\
              expand("data/{{normalized_dataset,.*?[a-z]}}/csv/delaunay/delaunay_{algorithm}.csv",
                     algorithm=["average_distance", "total_distance", "cartesian_product",
                                "number_instances", "frequency_instances"])
    params:
         snakefile="nodes/meta_workflows/structure_based_encodings/structure_based_encodings_windowed.smk",
         configfile="nodes/meta_workflows/structure_based_encodings/config.yaml"
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

#######################################################################################################################
############################################### COLLECT ENCODINGS #####################################################
#######################################################################################################################

rule utils_collect_encodings_bachem:
    input:
         sequence_based_encodings_in=\
             expand(rules.meta_workflow_sequence_based_encodings.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if "complete" not in ds and "bachem" in ds]),
         structure_based_encodings_in=\
             expand(rules.meta_workflow_structure_based_encodings.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if "complete" not in ds and int(ds.split("_")[-1]) >= 12 and "bachem" in ds]) + \
             expand(rules.meta_workflow_structure_based_encodings_windowed.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if "complete" in ds and "bachem" in ds])
    output:
         sequence_based_encodings_out=\
             directory(f"data/bachem/csv/sequence_based/"),
         structure_based_encodings_out=\
             directory(f"data/bachem/csv/structure_based/"),
         csv_dir_out=\
             directory(f"data/bachem/csv/non_empty/all/")
    params:
         snakefile="nodes/utils/collect_encodings/Snakefile",
         configfile="nodes/utils/collect_encodings/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=DATASETS, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_collect_encodings_protease:
    input:
         sequence_based_encodings_in=\
             expand(rules.meta_workflow_sequence_based_encodings.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if "complete" not in ds and "protease" in ds]),
         structure_based_encodings_in=\
             expand(rules.meta_workflow_structure_based_encodings.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if "complete" not in ds and int(ds.split("_")[-1]) >= 12 and "protease" in ds]) + \
             expand(rules.meta_workflow_structure_based_encodings_windowed.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if "complete" in ds and "protease" in ds])
    output:
         sequence_based_encodings_out=\
             directory(f"data/protease/csv/sequence_based/"),
         structure_based_encodings_out=\
             directory(f"data/protease/csv/structure_based/"),
         csv_dir_out=\
             directory(f"data/protease/csv/non_empty/all/")
    params:
         snakefile="nodes/utils/collect_encodings/Snakefile",
         configfile="nodes/utils/collect_encodings/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=DATASETS, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule plot_empty_datasets_protease:
    input:
         sequence_based_encodings=f"data/protease/csv/sequence_based/",
         structure_based_encodings=f"data/protease/csv/structure_based/"
    output:
         png_out=f"data/protease/plots/filtered_datasets.png"
    params:
         snakefile="nodes/plots/empty_datasets/Snakefile",
         configfile="nodes/plots/empty_datasets/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=[ds for ds in DATASETS if "protease" in ds], cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule plot_empty_datasets_bachem:
    input:
         sequence_based_encodings=f"data/bachem/csv/sequence_based/",
         structure_based_encodings=f"data/bachem/csv/structure_based/"
    output:
         png_out=f"data/bachem/plots/filtered_datasets.png"
    params:
         snakefile="nodes/plots/empty_datasets/Snakefile",
         configfile="nodes/plots/empty_datasets/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=[ds for ds in DATASETS if "bachem" in ds], cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

########################################################################################################################
################################################ MACHINE LEARNING ######################################################
########################################################################################################################

# rule machine_learning_hold_out_datasets_bachem:
#     input:
#          csv_dir_in=f"data/bachem/csv/non_empty/all/"
#     output:
#          csv_train_dir_out=directory(f"data/bachem/csv/non_empty_train/"),
#          csv_val_1_dir_out=directory(f"data/bachem/csv/non_empty_val_1/"),
#          csv_val_2_dir_out=directory(f"data/bachem/csv/non_empty_val_2/"),
#          csv_test_dir_out=directory(f"data/bachem/csv/non_empty_test/")
#     params:
#          snakefile="nodes/machine_learning/hold_out_datasets/Snakefile",
#          configfile="nodes/machine_learning/hold_out_datasets/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
#
# rule machine_learning_hold_out_datasets_protease:
#     input:
#          csv_dir_in=f"data/protease/csv/non_empty/all/"
#     output:
#          csv_train_dir_out=directory(f"data/protease/csv/non_empty_train/"),
#          csv_val_1_dir_out=directory(f"data/protease/csv/non_empty_val_1/"),
#          csv_val_2_dir_out=directory(f"data/protease/csv/non_empty_val_2/"),
#          csv_test_dir_out=directory(f"data/protease/csv/non_empty_test/")
#     params:
#          snakefile="nodes/machine_learning/hold_out_datasets/Snakefile",
#          configfile="nodes/machine_learning/hold_out_datasets/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
#
# rule machine_learning_top_encodings_bachem:
#     input:
#          train_dir_in=f"data/bachem/csv/non_empty_train/",
#          val_dir_in=f"data/bachem/csv/non_empty_val_1/"
#     output:
#          top_encodings_out="data/bachem/machine_learning/top_encodings.csv",
#          phi_correlation_out=f"data/bachem/machine_learning/phi_correlation.csv",
#     params:
#          snakefile="nodes/machine_learning/top_encodings/Snakefile",
#          configfile="nodes/machine_learning/top_encodings/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
#
# rule machine_learning_top_encodings_protease:
#     input:
#          train_dir_in=f"data/protease/csv/non_empty_train/",
#          val_dir_in=f"data/protease/csv/non_empty_val_1/"
#     output:
#          top_encodings_out="data/protease/machine_learning/top_encodings.csv",
#          phi_correlation_out=f"data/protease/machine_learning/phi_correlation.csv"
#     params:
#          snakefile="nodes/machine_learning/top_encodings/Snakefile",
#          configfile="nodes/machine_learning/top_encodings/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
#
# # rule machine_learning_best_ensemble_bachem:
# #     input:
# #          train_dirs_in=[f"data/bachem/csv/non_empty_train/", f"data/bachem/csv/non_empty_val_1/"],
# #          val_dir_in=f"data/bachem/csv/non_empty_val_2/",
# #          test_dir_in=f"data/bachem/csv/non_empty_train/",
# #          top_encodings_in="data/bachem/machine_learning/top_encodings.csv"
# #     output:
# #          ensemble_out=f"data/bachem/machine_learning/diff_coupl_clf.jl"
# #     run:
# #          # TODO https://scikit-learn.org/stable/modules/ensemble.html#voting-classifier
# #          pass
#
# rule machine_learning_best_ensemble_protease:
#     input:
#          train_dirs_in=[f"data/protease/csv/non_empty_train/", f"data/protease/csv/non_empty_val_1/"],
#          val_dir_in=f"data/protease/csv/non_empty_val_2/",
#          test_dir_in=f"data/protease/csv/non_empty_test/",
#          phi_correlation_in=f"data/protease/machine_learning/phi_correlation.csv"
#     output:
#          ensemble_validation_out=\
#              "data/protease/machine_learning/ensemble_results_validated.csv",
#          ensemble_validation_tuned_hp_out=\
#              "data/protease/machine_learning/ensemble_results_validated_tuned_hp.csv"
#     params:
#          snakefile="nodes/machine_learning/best_ensemble/Snakefile",
#          configfile="nodes/machine_learning/best_ensemble/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")
#
#
#
#
