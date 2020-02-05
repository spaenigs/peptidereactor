from proteinreactor.workflow_executer import WorkflowExecuter, MetaWorkflowExecuter
import pandas as pd

WINDOW_LENGTHS = [8,12]#[8,11,15,20]
DATASETS = \
    ["protease_window_length_8", "protease_window_length_8_complete"] #+ \
    # expand(["bachem_window_length_{window_length}", "bachem_window_length_{window_length}_complete"],
    #        window_length=WINDOW_LENGTHS)

CORES = int(config["cores"])

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

rule all:
    input:
         # "data/protease_window_length_8/seqs.fasta",
         # "data/protease_window_length_8/classes.txt",
         # "data/protease_window_length_8_complete/seqs_orig.fasta",
         # "data/protease_window_length_8_complete/classes_orig.txt",
         # "data/protease_window_length_8_complete/seqs.fasta",
         # "data/protease_window_length_8_complete/classes.yaml",
         # "data/protease_window_length_8_complete/classes.txt",
         # expand("data/{normalized_dataset}/csv/zscale.csv", normalized_dataset=DATASETS[0]),
         # expand("data/{normalized_dataset}/csv/tpc.csv", normalized_dataset=DATASETS[0]),
         # expand("data/{normalized_dataset}/csv/asa.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/ta.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/ssec.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/sseb.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/disorder.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/disorderb.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/disorderc.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/qsar.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
         #        normalized_dataset=DATASETS[1], distance=[0,3,6,9,12]),
         # expand("data/{normalized_dataset}/csv/distance_distribution.csv", normalized_dataset=DATASETS[1]),
         # expand("data/{normalized_dataset}/csv/delaunay/delaunay_{algorithm}.csv",
         #        normalized_dataset=DATASETS[1],
         #        algorithm=["average_distance", "total_distance", "cartesian_product",
         #                   "number_instances", "frequency_instances"])
         f"data/bachem/plots/filtered_datasets.png",
         f"data/bachem/machine_learning/top_encodings.csv"

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
         # gtpc_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/gtpc.csv",
         # gdpc_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/gdpc.csv",
         # gaac_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/gaac.csv",
         # dpc_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/dpc.csv",
         # dde_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/dde.csv",
         # ctdt_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ctdt.csv",
         # ctdd_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ctdd.csv",
         # ctdc_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ctdc.csv",
         # blosum62_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/blosum62.csv",
         # binary_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/binary.csv",
         # aac_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/aac.csv",
         # ctriad_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ctriad.csv",
         # blomap_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/blomap.csv",
         # egaac_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/egaac/egaac_window_{window_val}.csv", window_val=list(range(1, 31))),
         # aaindex_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/aaindex/aaindex_{aaindex}.csv", aaindex=get_aaindex()),
         # fft_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/fft/fft_{aaindex}.csv", aaindex=get_aaindex()),
         # waac_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/waac/waac_{aaindex}.csv", aaindex=get_aaindex()),
         # flgc_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/flgc/flgc_{aaindex}.csv", aaindex=get_aaindex()),
         # fldpc_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/fldpc/fldpc_{aaindex}.csv", aaindex=get_aaindex()),
         # ngram_a2_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ngram_a2/ngram_a2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         # ngram_a3_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ngram_a3/ngram_a3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         # ngram_e2_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ngram_e2/ngram_e2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         # ngram_e3_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ngram_e3/ngram_e3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         # ngram_s2_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ngram_s2/ngram_s2_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         # ngram_s3_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ngram_s3/ngram_s3_{dim}.csv", dim=[1, 5, 20, 50, 100, 200, 300]),
         # cgr_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
         #             resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),
         # distance_frequency_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/distance_frequency/distance_frequency_dn_{nterminal}_dc_{cterminal}.csv",
         #             nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),
         # cksaagp_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/cksaagp/cksaagp_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         # socnumber_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/socnumber/socnumber_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         # qsorder_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/qsorder/qsorder_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         # nmbroto_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         # moran_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/moran/moran_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         # ksctriad_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/ksctriad/ksctriad_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         # geary_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/geary/geary_nlag_{nlag_val}.csv", nlag_val=list(range(1, 31))),
         # eaac_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/eaac/eaac_window_{window_val}.csv", window_val=list(range(1, 31))),
         # cksaap=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/cksaap/cksaap_gap_{gap_val}.csv", gap_val=list(range(1, 31))),
         # apaac_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/apaac/apaac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 31))),
         # paac_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/paac/paac_lambda_{lambda_val}.csv", lambda_val=list(range(1, 31))),
         # psekraac_type16_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type16/psekraac_type16_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type15_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type15/psekraac_type15_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type14_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type14/psekraac_type14_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type13_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type13/psekraac_type13_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type12_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type12/psekraac_type12_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type11_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type11/psekraac_type11_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type10_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type10/psekraac_type10_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type9_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type9/psekraac_type9_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type8_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type8/psekraac_type8_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type7_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type7/psekraac_type7_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type6C_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type6C/psekraac_type6C_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type6B_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type6B/psekraac_type6B_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type6A_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type6A/psekraac_type6A_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type5_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type5/psekraac_type5_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type4_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type4/psekraac_type4_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type3B_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type3B/psekraac_type3B_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type3A_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type3A/psekraac_type3A_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type2_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type2/psekraac_type2_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # psekraac_type1_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/psekraac_type1/psekraac_type1_"
         #             "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #             sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #             ktuple_val=list(range(1,4)), lambda_val=list(range(1,4)))
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
         fasta_anno_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_seqs.fasta",
         classes_anno="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?(1[2-9]|2\d)}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?(1[2-9]|2\d)}/pdb/"),
         # pssm_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/pssm.csv",
         asa_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/asa.csv",
         # ta_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ta.csv",
         # ssec_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/ssec.csv",
         # sseb_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/sseb.csv",
         # disorder_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/disorder.csv",
         # disorderb_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/disorderb.csv",
         # disorderc_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/disorderc.csv",
         # qsar_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/qsar.csv",
         # electrostatic_hull_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
         #             distance=[0,3,6,9,12]),
         # distance_distribution_out="data/{normalized_dataset,.*?(1[2-9]|2\d)}/csv/distance_distribution.csv",
         # delaunay_out=\
         #      expand("data/{{normalized_dataset,.*?(1[2-9]|2\d)}}/csv/delaunay/delaunay_{algorithm}.csv",
         #             algorithm=["average_distance", "total_distance", "cartesian_product",
         #                        "number_instances", "frequency_instances"])
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
         classes_in="data/{normalized_dataset}/classes.yaml",
         # fasta_anno_out="data/{normalized_dataset}/annotated_seqs.fasta",
         # classes_anno_idx_out="data/{normalized_dataset}/annotated_classes.txt",
         # fasta_anno_msa_out="data/{normalized_dataset}/annotated_seqs_msa.fasta",
         # profile_dir=directory("data/{normalized_dataset}/profile/"),
         # fasta_anno_pdbs_out="data/{normalized_dataset}/annotated_pdbs_seqs.fasta",
         # classes_anno_pdbs_idx_out="data/{normalized_dataset}/annotated_pdbs_classes.txt",
         # pdb_out=directory("data/{normalized_dataset}/pdb/"),
    output:
         fasta_anno_out="data/{normalized_dataset,.*?[a-z]}/annotated_seqs.fasta",
         classes_anno_idx_out="data/{normalized_dataset,.*?[a-z]}/annotated_classes.txt",
         fasta_anno_msa_out="data/{normalized_dataset,.*?[a-z]}/annotated_seqs_msa.fasta",
         profile_dir=directory("data/{normalized_dataset,.*?[a-z]}/profile/"),
         fasta_anno_pdbs_out="data/{normalized_dataset,.*?[a-z]}/annotated_pdbs_seqs.fasta",
         classes_anno_pdbs_idx_out="data/{normalized_dataset,.*?[a-z]}/annotated_pdbs_classes.txt",
         pdb_out=directory("data/{normalized_dataset,.*?[a-z]}/pdb/"),
         # asa_out="data/{normalized_dataset,.*?[a-z]}/csv/asa.csv",
         # ta_out="data/{normalized_dataset,.*?[a-z]}/csv/ta.csv",
         # ssec_out="data/{normalized_dataset,.*?[a-z]}/csv/ssec.csv",
         # sseb_out="data/{normalized_dataset,.*?[a-z]}/csv/sseb.csv",
         # disorder_out="data/{normalized_dataset,.*?[a-z]}/csv/disorder.csv",
         # disorderb_out="data/{normalized_dataset,.*?[a-z]}/csv/disorderb.csv",
         # disorderc_out="data/{normalized_dataset,.*?[a-z]}/csv/disorderc.csv",
         # qsar_out="data/{normalized_dataset,.*?[a-z]}/csv/qsar.csv",
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

########################################################################################################################
################################################ COLLECT ENCODINGS #####################################################
########################################################################################################################

rule collect_encodings:
    input:
         sequence_based_encodings=\
             expand(rules.meta_workflow_sequence_based_encodings.output,
                    normalized_dataset=[ds for ds in DATASETS
                                        if ds.split("_")[-1] != "complete"]),
         structure_based_encodings=\
             expand(rules.meta_workflow_structure_based_encodings.output[7:],
                    normalized_dataset=[ds for ds in DATASETS
                                        if ds.split("_")[-1] != "complete" and int(ds.split("_")[-1]) >= 12]) + \
             expand(rules.meta_workflow_structure_based_encodings_windowed.output[7:],
                    normalized_dataset=[ds for ds in DATASETS
                                        if ds.split("_")[-1] == "complete"])
    output:
         sequence_based_encodings=\
             directory(f"data/temp/datasets/csv/sequence_based/"),
         structure_based_encodings=\
             directory(f"data/temp/datasets/csv/structure_based/")
    run:
         import os
         def copy(from_files, target_dir):
            for i in from_files:
                dataset = i.split("/")[1]
                shell(f"cp {i} {target_dir}{dataset}_{os.path.basename(i)}")
         copy(list(input.sequence_based_encodings), str(output.sequence_based_encodings))
         copy(list(input.structure_based_encodings), str(output.structure_based_encodings))

rule remove_empty_datasets:
    input:
         f"data/temp/datasets/csv/sequence_based/",
         f"data/temp/datasets/csv/structure_based/"
    output:
         directory(f"data/temp/datasets/csv_non_empty/")
    run:
         from glob import glob

         for csv_path in glob(input[0] + "*.csv") + glob(input[1] + "*.csv"):
             if os.path.getsize(csv_path) == 0:
                 continue
             else:
                 shell(f"cp {csv_path} {str(output)}")

rule plot_empty_datasets:
    input:
         sequence_based_encodings=f"data/temp/datasets/csv/sequence_based/",
         structure_based_encodings=f"data/temp/datasets/csv/structure_based/"
    output:
         png_out=f"data/bachem/plots/filtered_datasets.png"
    params:
         snakefile="nodes/plots/empty_datasets/Snakefile",
         configfile="nodes/plots/empty_datasets/config.yaml"
    priority:
        50
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, datasets=DATASETS) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

########################################################################################################################
################################################ MACHINE LEARNING ######################################################
########################################################################################################################

# TODO keep final test set
rule hold_out_data:
    input:
         f"data/temp/datasets/csv_non_empty/"
    output:
         directory(f"data/temp/datasets/csv_non_empty_train/"),
         directory(f"data/temp/datasets/csv_non_empty_test/")
    run:
         from sklearn.model_selection import train_test_split
         from glob import glob

         for csv_path in glob(str(input) + "*.csv"):
             csv_name = os.path.basename(csv_path)
             df = pd.read_csv(csv_path, index_col=0)
             X_train, X_test, y_train, y_test = \
                 train_test_split(df.iloc[:,:-1], df["y"], test_size=0.1, random_state=42, stratify=df["y"])
             X_train["y"], X_test["y"] = y_train, y_test
             X_train.to_csv(str(output[0]) + csv_name)
             X_test.to_csv(str(output[1]) + csv_name)

rule machine_learning_top_encodings:
    input:
         csv_dir_in=f"data/temp/datasets/csv_non_empty_train/"
    output:
         csv_out=f"data/bachem/machine_learning/top_encodings.csv"
    params:
         snakefile="nodes/machine_learning/top_encodings/Snakefile",
         configfile="nodes/machine_learning/top_encodings/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# TODO https://scikit-learn.org/stable/modules/ensemble.html#voting-classifier

