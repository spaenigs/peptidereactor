from proteinreactor.workflow_executer import WorkflowExecuter, MetaWorkflowExecuter
import pandas as pd

DATASET = "bachem_test"
# DATASETS = expand("bachem_window_length_{window_length}", window_length=[8,11,15,20]) + ["protease"]
DATASETS = ["bachem_test_window_length_8", "bachem_test_window_length_8_complete"]
CORES = int(config["cores"])

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

rule all:
    input:
         # expand(f"data/{DATASET}_window_length_{{window_length}}_complete/seqs.fasta", window_length=[8,11,15,20]),
         # f"data/{DATASET}_window_length_7_complete/seqs.fasta",
         # f"data/{DATASET}_window_length_7_complete/classes.yaml"
         # expand(f"data/{DATASET}_window_length_{{window_length}}_complete/classes.yaml", window_length=[8,11,15,20])
         f"data/{DATASET}/plots/filtered_datasets.png",
         f"data/{DATASET}/machine_learning/top_encodings.csv"

########################################################################################################################
############################################## DATASET CREATION ########################################################
########################################################################################################################

rule utils_sliding_windows:
    input:
         series_in=f"data/{DATASET}/series.yaml"
    output:
         fastas_out=f"data/{DATASET}_window_length_8/seqs.fasta",
         classes_out=f"data/{DATASET}_window_length_8/classes.txt",
         # fastas_out=expand(f"data/{DATASET}_window_length_{{window_length}}/seqs.fasta", window_length=[0,8,11,15,20]),
         # classes_out=expand(f"data/{DATASET}_window_length_{{window_length}}/classes.txt", window_length=[0,8,11,15,20])
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset=DATASET) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_sliding_windows_complete:
    input:
         series_in=f"data/{DATASET}/series.yaml"
    output:
         fastas_out=f"data/{DATASET}_window_length_8_complete/seqs.fasta",
         classes_out=f"data/{DATASET}_window_length_8_complete/classes.yaml",
         classes_idx_out=f"data/{DATASET}_window_length_8_complete/classes.txt"
         # fastas_out=expand(f"data/{DATASET}_window_length_{{window_length}}_complete/seqs.fasta", window_length=[8,11,15,20]),
         # classes_out=expand(f"data/{DATASET}_window_length_{{window_length}}_complete/classes.yaml", window_length=[8,11,15,20]),
         # classes_idx_out=expand(f"data/{DATASET}_window_length_{{window_length}}_complete/classes.txt", window_length=[8,11,15,20])
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows_complete.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset=DATASET) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

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
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset=DATASET) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

########################################################################################################################
############################################ SEQUENCE-BASED ENCODINGS #################################################
########################################################################################################################

rule meta_workflow_sequence_based_encodings:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
         fasta_anno_in="data/{normalized_dataset}/annotated_seqs.fasta",
         classes_anno_in="data/{normalized_dataset}/annotated_classes.txt",
         profile_dir="data/{normalized_dataset}/profile/"
    output:
         pssm_out="data/{normalized_dataset}/csv/pssm.csv",
         zscale_out="data/{normalized_dataset}/csv/zscale.csv",
         tpc_out="data/{normalized_dataset}/csv/tpc.csv",
         gtpc_out="data/{normalized_dataset}/csv/gtpc.csv",
         gdpc_out="data/{normalized_dataset}/csv/gdpc.csv",
         gaac_out="data/{normalized_dataset}/csv/gaac.csv",
         dpc_out="data/{normalized_dataset}/csv/dpc.csv",
         dde_out="data/{normalized_dataset}/csv/dde.csv",
         ctdt_out="data/{normalized_dataset}/csv/ctdt.csv",
         ctdd_out="data/{normalized_dataset}/csv/ctdd.csv",
         ctdc_out="data/{normalized_dataset}/csv/ctdc.csv",
         blosum62_out="data/{normalized_dataset}/csv/blosum62.csv",
         binary_out="data/{normalized_dataset}/csv/binary.csv",
         aac_out="data/{normalized_dataset}/csv/aac.csv",
         ctriad_out="data/{normalized_dataset}/csv/ctriad.csv",
         blomap_out="data/{normalized_dataset}/csv/blomap.csv",
         egaac_out=expand("data/{{normalized_dataset}}/csv/egaac/egaac_window_{window_val}.csv", window_val=list(range(1, 31))),
         aaindex_out=expand("data/{{normalized_dataset}}/csv/aaindex/aaindex_{aaindex}.csv", aaindex=get_aaindex()),
         fft_out=expand("data/{{normalized_dataset}}/csv/fft/fft_{aaindex}.csv", aaindex=get_aaindex()),
         waac_out=expand("data/{{normalized_dataset}}/csv/waac/waac_{aaindex}.csv", aaindex=get_aaindex()),
         flgc_out=expand("data/{{normalized_dataset}}/csv/flgc/flgc_{aaindex}.csv", aaindex=get_aaindex()),
         fldpc_out=expand("data/{{normalized_dataset}}/csv/fldpc/fldpc_{aaindex}.csv", aaindex=get_aaindex()),
         ngram_a2_out=expand("data/{{normalized_dataset}}/csv/ngram_a2/ngram_a2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_a3_out=expand("data/{{normalized_dataset}}/csv/ngram_a3/ngram_a3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e2_out=expand("data/{{normalized_dataset}}/csv/ngram_e2/ngram_e2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_e3_out=expand("data/{{normalized_dataset}}/csv/ngram_e3/ngram_e3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s2_out=expand("data/{{normalized_dataset}}/csv/ngram_s2/ngram_s2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         ngram_s3_out=expand("data/{{normalized_dataset}}/csv/ngram_s3/ngram_s3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         cgr_out=expand("data/{{normalized_dataset}}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv", resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),
         distance_frequency_out=expand("data/{{normalized_dataset}}/csv/distance_frequency/distance_frequency_dn_{nterminal}_dc_{cterminal}.csv", nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),
         cksaagp_out=expand("data/{{normalized_dataset}}/csv/cksaagp/cksaagp_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         socnumber_out=expand("data/{{normalized_dataset}}/csv/socnumber/socnumber_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         qsorder_out=expand("data/{{normalized_dataset}}/csv/qsorder/qsorder_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         nmbroto_out=expand("data/{{normalized_dataset}}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         moran_out=expand("data/{{normalized_dataset}}/csv/moran/moran_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         ksctriad_out=expand("data/{{normalized_dataset}}/csv/ksctriad/ksctriad_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         geary_out=expand("data/{{normalized_dataset}}/csv/geary/geary_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         eaac_out=expand("data/{{normalized_dataset}}/csv/eaac/eaac_window_{window_val}.csv", window_val=list(range(1, 31))),
         cksaap=expand("data/{{normalized_dataset}}/csv/cksaap/cksaap_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         apaac_out=expand("data/{{normalized_dataset}}/csv/apaac/apaac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 31))),
         paac_out=expand("data/{{normalized_dataset}}/csv/paac/paac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 31))),
         psekraac_type16_out=expand("data/{{normalized_dataset}}/csv/psekraac_type16/psekraac_type16_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type15_out=expand("data/{{normalized_dataset}}/csv/psekraac_type15/psekraac_type15_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type14_out=expand("data/{{normalized_dataset}}/csv/psekraac_type14/psekraac_type14_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type13_out=expand("data/{{normalized_dataset}}/csv/psekraac_type13/psekraac_type13_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type12_out=expand("data/{{normalized_dataset}}/csv/psekraac_type12/psekraac_type12_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type11_out=expand("data/{{normalized_dataset}}/csv/psekraac_type11/psekraac_type11_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type10_out=expand("data/{{normalized_dataset}}/csv/psekraac_type10/psekraac_type10_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type9_out=expand("data/{{normalized_dataset}}/csv/psekraac_type9/psekraac_type9_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type8_out=expand("data/{{normalized_dataset}}/csv/psekraac_type8/psekraac_type8_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type7_out=expand("data/{{normalized_dataset}}/csv/psekraac_type7/psekraac_type7_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type6C_out=expand("data/{{normalized_dataset}}/csv/psekraac_type6C/psekraac_type6C_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[5], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type6B_out=expand("data/{{normalized_dataset}}/csv/psekraac_type6B/psekraac_type6B_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[5], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type6A_out=expand("data/{{normalized_dataset}}/csv/psekraac_type6A/psekraac_type6A_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type5_out=expand("data/{{normalized_dataset}}/csv/psekraac_type5/psekraac_type5_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type4_out=expand("data/{{normalized_dataset}}/csv/psekraac_type4/psekraac_type4_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type3B_out=expand("data/{{normalized_dataset}}/csv/psekraac_type3B/psekraac_type3B_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type3A_out=expand("data/{{normalized_dataset}}/csv/psekraac_type3A/psekraac_type3A_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type2_out=expand("data/{{normalized_dataset}}/csv/psekraac_type2/psekraac_type2_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20], ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         psekraac_type1_out=expand("data/{{normalized_dataset}}/csv/psekraac_type1/psekraac_type1_subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv", sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)), ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
    params:
         snakefile="nodes/meta_workflows/sequence_based_encodings/Snakefile",
         configfile="nodes/meta_workflows/sequence_based_encodings/config.yaml"
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

########################################################################################################################
############################################ STRUCTURE-BASED ENCODINGS #################################################
########################################################################################################################

rule meta_workflow_structure_based_encodings:
    input:
         fasta_in="data/{normalized_dataset}/seqs.fasta",
         fasta_msa_in="data/{normalized_dataset}/seqs_msa.fasta",
         classes_in="data/{normalized_dataset}/classes.txt",
    output:
         fasta_anno_out="data/{normalized_dataset,.*?\d}/annotated_seqs.fasta",
         classes_anno="data/{normalized_dataset,.*?\d}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?\d}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?\d}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?\d}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_out="data/{normalized_dataset,.*?\d}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?\d}/pdb/"),
         asa_out="data/{normalized_dataset,.*?\d}/csv/asa.csv",
         ta_out="data/{normalized_dataset,.*?\d}/csv/ta.csv",
         ssec_out="data/{normalized_dataset,.*?\d}/csv/ssec.csv",
         sseb_out="data/{normalized_dataset,.*?\d}/csv/sseb.csv",
         disorder_out="data/{normalized_dataset,.*?\d}/csv/disorder.csv",
         disorderb_out="data/{normalized_dataset,.*?\d}/csv/disorderb.csv",
         disorderc_out="data/{normalized_dataset,.*?\d}/csv/disorderc.csv",
         qsar_out="data/{normalized_dataset,.*?\d}/csv/qsar.csv",
         electrostatic_hull_out=\
              expand("data/{{normalized_dataset,.*?\d}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                     distance=[0,3,6,9,12]),
         distance_distribution_out="data/{normalized_dataset,.*?\d}/csv/distance_distribution.csv",
         delaunay_out=\
              expand("data/{{normalized_dataset,.*?\d}}/csv/delaunay/delaunay_{algorithm}.csv",
                     algorithm=["average_distance", "total_distance", "cartesian_product",
                                "number_instances", "frequency_instances"])
    params:
         snakefile="nodes/meta_workflows/structure_based_encodings/structure_based_encodings.smk",
         configfile="nodes/meta_workflows/structure_based_encodings/config.yaml"
    run:
         with MetaWorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule meta_workflow_structure_based_encodings_windowed:
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
         pdb_out=directory("data/{normalized_dataset,.*?[a-z]}/pdb/"),
         asa_out="data/{normalized_dataset,.*?[a-z]}/csv/asa.csv",
         ta_out="data/{normalized_dataset,.*?[a-z]}/csv/ta.csv",
         ssec_out="data/{normalized_dataset,.*?[a-z]}/csv/ssec.csv",
         sseb_out="data/{normalized_dataset,.*?[a-z]}/csv/sseb.csv",
         disorder_out="data/{normalized_dataset,.*?[a-z]}/csv/disorder.csv",
         disorderb_out="data/{normalized_dataset,.*?[a-z]}/csv/disorderb.csv",
         disorderc_out="data/{normalized_dataset,.*?[a-z]}/csv/disorderc.csv",
         qsar_out="data/{normalized_dataset,.*?[a-z]}/csv/qsar.csv",
         electrostatic_hull_out=\
              expand("data/{{normalized_dataset,.*?[a-z]$}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
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

########################################################################################################################
################################################ MACHINE LEARNING ######################################################
########################################################################################################################

rule collect_encodings:
    input:
         sequence_based_encodings=\
             expand("data/{normalized_dataset}/csv/aaindex/aaindex_{aaindex}.csv",
                    normalized_dataset=filter(lambda ds: "complete" not in ds, DATASETS), aaindex=get_aaindex()) +
             # expand("data/{normalized_dataset}/csv/apaac/apaac_lambda_{lambda_val}.csv",
             #        normalized_dataset=DATASETS, lambda_val=list(range(1, 31))) +
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
             # expand("data/{normalized_dataset}/csv/aac.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/binary.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/blosum62.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/ctdc.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/ctdd.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/ctdt.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/ctriad.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/dde.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/dpc.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/gaac.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/gdpc.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/gtpc.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/tpc.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/zscale.csv", normalized_dataset=DATASETS) +
             expand("data/{normalized_dataset}/csv/blomap.csv",
                    normalized_dataset=filter(lambda ds: "complete" not in ds, DATASETS)),
         structure_based_encodings=\
             expand("data/{normalized_dataset}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                    normalized_dataset=DATASETS, distance=[0,3,6,9,12]) +
             # expand("data/{normalized_dataset}/csv/asa.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/delaunay/delaunay_{algorithm}.csv",
             #        normalized_dataset=DATASETS,
             #        algorithm=["average_distance", "total_distance", "cartesian_product",
             #                   "number_instances", "frequency_instances"]) +
             # expand("data/{normalized_dataset}/csv/disorder.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/disorderb.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/disorderc.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/pssm.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/sseb.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/ssec.csv", normalized_dataset=DATASETS) +
             # expand("data/{normalized_dataset}/csv/ta.csv", normalized_dataset=DATASETS) +
             expand("data/{normalized_dataset}/csv/distance_distribution.csv", normalized_dataset=DATASETS) +
             expand("data/{normalized_dataset}/csv/qsar.csv", normalized_dataset=DATASETS) +
             expand("data/{normalized_dataset}/csv/asa.csv", normalized_dataset=DATASETS)
    output:
         # directory("data/temp/{normalized_dataset}/")
         sequence_based_encodings=\
             temp(expand("data/temp/{normalized_dataset}/csv/aaindex/aaindex_{aaindex}.csv",
                    normalized_dataset=filter(lambda ds: "complete" not in ds, DATASETS), aaindex=get_aaindex())) +
             # temp(expand("data/temp/{normalized_dataset}/csv/apaac/apaac_lambda_{lambda_val}.csv",
             #        normalized_dataset=DATASETS, lambda_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/cksaagp/cksaagp_gap_{gap_val}.csv",
             #        normalized_dataset=DATASETS, gap_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/socnumber/socnumber_nlag_{nlag_val}.csv",
             #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/qsorder/qsorder_nlag_{nlag_val}.csv",
             #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv",
             #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/moran/moran_nlag_{nlag_val}.csv",
             #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ksctriad/ksctriad_gap_{gap_val}.csv",
             #        normalized_dataset=DATASETS, gap_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/geary/geary_nlag_{nlag_val}.csv",
             #        normalized_dataset=DATASETS, nlag_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/eaac/eaac_window_{window_val}.csv",
             #        normalized_dataset=DATASETS, window_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/cksaap/cksaap_gap_{gap_val}.csv",
             #        normalized_dataset=DATASETS, gap_val=list(range(1, 31)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/paac/paac_lambda_{lambda_val}.csv",
             #        normalized_dataset=DATASETS, lambda_val=list(range(1, 31)))) +
             #
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type1/psekraac_type1_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type2/psekraac_type2_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type3A/psekraac_type3A_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type3B/psekraac_type3B_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type4/psekraac_type4_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type5/psekraac_type5_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type6A/psekraac_type6A_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type6B/psekraac_type6B_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type6C/psekraac_type6C_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type7/psekraac_type7_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type8/psekraac_type8_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type9/psekraac_type9_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type10/psekraac_type10_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type11/psekraac_type11_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type12/psekraac_type12_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type13/psekraac_type13_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type14/psekraac_type14_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type15/psekraac_type15_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             # temp(expand("data/temp/{normalized_dataset}/csv/psekraac_type16/psekraac_type16_"
             #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
             #        normalized_dataset=DATASETS,
             #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
             #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))) +
             #
             # temp(expand("data/temp/{normalized_dataset}/csv/fft/fft_{aaindex}.csv",
             #        normalized_dataset=DATASETS, aaindex=get_aaindex())) +
             #
             # temp(expand("data/temp/{normalized_dataset}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
             #        normalized_dataset=DATASETS,
             #        resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/waac/waac_{aaindex}.csv",
             #        normalized_dataset=DATASETS,
             #        aaindex=get_aaindex())) +
             # temp(expand("data/temp/{normalized_dataset}/csv/flgc/flgc_{aaindex}.csv",
             #        normalized_dataset=DATASETS,
             #        aaindex=get_aaindex())) +
             # temp(expand("data/temp/{normalized_dataset}/csv/fldpc/fldpc_{aaindex}.csv",
             #        normalized_dataset=DATASETS,
             #        aaindex=get_aaindex())) +
             #
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_a2/ngram_a2_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_a2/ngram_a2_lsv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_a2/ngram_a2_sv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_a3/ngram_a3_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_a3/ngram_a3_lsv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_a3/ngram_a3_sv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_e2/ngram_e2_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_e2/ngram_e2_lsv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_e2/ngram_e2_sv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_e3/ngram_e3_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_e3/ngram_e3_lsv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_e3/ngram_e3_sv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_s2/ngram_s2_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_s2/ngram_s2_lsv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_s2/ngram_s2_sv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_s3/ngram_s3_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_s3/ngram_s3_lsv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ngram_s3/ngram_s3_sv_{dim}.csv",
             #        normalized_dataset=DATASETS, dim=[1, 5, 20, 50, 100, 200, 300])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/aac.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/binary.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/blosum62.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ctdc.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ctdd.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ctdt.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ctriad.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/dde.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/dpc.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/gaac.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/gdpc.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/gtpc.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/tpc.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/zscale.csv", normalized_dataset=DATASETS)) +
             temp(expand("data/temp/{normalized_dataset}/csv/blomap.csv",
                         normalized_dataset=filter(lambda ds: "complete" not in ds, DATASETS))),
         structure_based_encodings=\
             temp(expand("data/temp/{normalized_dataset}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
                    normalized_dataset=DATASETS, distance=[0,3,6,9,12])) +
             temp(expand("data/temp/{normalized_dataset}/csv/asa.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/delaunay/delaunay_{algorithm}.csv",
             #        normalized_dataset=DATASETS,
             #        algorithm=["average_distance", "total_distance", "cartesian_product",
             #                   "number_instances", "frequency_instances"])) +
             # temp(expand("data/temp/{normalized_dataset}/csv/disorder.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/disorderb.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/disorderc.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/pssm.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/sseb.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ssec.csv", normalized_dataset=DATASETS)) +
             # temp(expand("data/temp/{normalized_dataset}/csv/ta.csv", normalized_dataset=DATASETS)) +
             temp(expand("data/temp/{normalized_dataset}/csv/distance_distribution.csv", normalized_dataset=DATASETS)) +
             temp(expand("data/temp/{normalized_dataset}/csv/qsar.csv", normalized_dataset=DATASETS))
    run:
         for i in list(input):
             target = i.replace("data", "data/temp")
             shell("cp {i} {target}")

rule plot_empty_datasets:
    input:
         sequence_based_encodings=rules.collect_encodings.output.sequence_based_encodings,
         structure_based_encodings=rules.collect_encodings.output.structure_based_encodings
    output:
         png_out=f"data/{DATASET}/plots/filtered_datasets.png"
    params:
         snakefile="nodes/plots/empty_datasets/Snakefile",
         configfile="nodes/plots/empty_datasets/config.yaml"
    priority:
        50
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=DATASETS) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule machine_learning_top_encodings:
    input:
         csvs_in=\
            rules.collect_encodings.output.sequence_based_encodings +
            rules.collect_encodings.output.structure_based_encodings
    output:
         csv_out=f"data/{DATASET}/machine_learning/top_encodings.csv"
    params:
         snakefile="nodes/machine_learning/top_encodings/Snakefile",
         configfile="nodes/machine_learning/top_encodings/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=DATASETS) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# https://scikit-learn.org/stable/modules/ensemble.html#voting-classifier