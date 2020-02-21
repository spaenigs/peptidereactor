from proteinreactor.workflow_executer import WorkflowExecuter, MetaWorkflowExecuter
from more_itertools import partition
import joblib
import re
import pandas as pd

CORES = config["cores"]

# bachem_profiles.smk
# bachem_encode.smk
# bachem_train.smk
# bachem_predict.smk

# [x] 1) take input sequence and split
# [x] 2) unpack model and get respective encodings
# [X] 3) if necessary, compute profile
# [x] 4) encoded splitted sequences
# [ ] 5) predict in model

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

def get_encoding_information():
    d = joblib.load("data/bachem/models/best_model_0.joblib")
    all, complete = partition(lambda ds: "complete" in ds, list(d.keys())[:3])
    tmp = {"datasets": [],
           "window_lengths": [],
           "encoding_names": [],
           "full_names": list(d.keys())[:3]}
    res = {"all": tmp, "complete": tmp}
    for dataset, encoding in [name.split("-") for name in all]:
        res["all"]["window_lengths"] += [int(re.search("(\d+)", dataset).group(0))]
        res["all"]["encoding_names"] += [encoding]
    for dataset, encoding in [name.split("-") for name in all]:
        res["complete"]["window_lengths"] += [int(re.search("(\d+)", dataset).group(0))]
        res["complete"]["encoding_names"] += [encoding]
    return res

ENCODING_PROFILE = get_encoding_information()

rule utils_sliding_windows:
    input:
         series_in=f"data/bachem/series_predict.yaml"
    output:
         fastas_out=expand(f"data/bachem_window_length_{{window_length}}/seqs.fasta",
                           window_length=ENCODING_PROFILE["all"]["window_lengths"]),
         classes_out=expand(f"data/bachem_window_length_{{window_length}}/classes.txt",
                            window_length=ENCODING_PROFILE["all"]["window_lengths"]),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule utils_sliding_windows_complete:
    input:
         series_in=f"data/bachem/series.yaml"
    output:
         fastas_out=\
             expand(f"data/bachem_window_length_{{window_length}}_complete/seqs.fasta",
                    window_length=ENCODING_PROFILE["complete"]["window_lengths"]),
         classes_out=\
              expand(f"data/bachem_window_length_{{window_length}}_complete/classes.yaml",
                     window_length=ENCODING_PROFILE["complete"]["window_lengths"]),
         classes_idx_out=\
              expand(f"data/bachem_window_length_{{window_length}}_complete/classes.txt",
                     window_length=ENCODING_PROFILE["complete"]["window_lengths"]),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows_complete.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
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

# compute profile only if necessary
encodings_all = \
    ["asa", "ta", "ssec", "sseb", "disorder", "disorderb", "disorderc",
     "qsar", "electrostatic_hull", "distance_distribution", "delaunay"]
filered_encodings_all = \
    [e for e in encodings_all if e in ENCODING_PROFILE["all"]["encoding_names"][0]] + \
    [e for e in encodings_all if e in ENCODING_PROFILE["all"]["encoding_names"][1]] + \
    [e for e in encodings_all if e in ENCODING_PROFILE["all"]["encoding_names"][2]]
if len(filered_encodings_all) > 0:
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

encodings_complete = \
    ["asa", "ta", "ssec", "sseb", "disorder", "disorderb", "disorderc",
     "qsar", "electrostatic_hull", "distance_distribution", "delaunay"]
filered_encodings_complete = \
    [e for e in encodings_complete if e in ENCODING_PROFILE["complete"]["encoding_names"][0]] + \
    [e for e in encodings_complete if e in ENCODING_PROFILE["complete"]["encoding_names"][1]] + \
    [e for e in encodings_complete if e in ENCODING_PROFILE["complete"]["encoding_names"][2]]
if len(filered_encodings_complete) > 0:
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

encodings_seq = ["zscale", "tpc", "gtpc", "gdpc", "gaac", "dpc", "dde", "ctdt", "ctdd", "ctdc", "blosum62", "binary",
                 "aac", "ctriad", "blomap", "egaac", "aaindex", "fft", "waac", "flgc", "fldpc", "ngram", "cgr",
                 "distance_frequency", "cksaap", "cksaagp", "socnumber", "qsorder", "nmbroto", "moran", "ksctriad",
                 "geary", "eaac", "apaac", "paac", "psekraac"]
filered_encodings_seq = \
    [e for e in encodings_seq if e in ENCODING_PROFILE["complete"]["encoding_names"][0]] + \
    [e for e in encodings_seq if e in ENCODING_PROFILE["complete"]["encoding_names"][1]] + \
    [e for e in encodings_seq if e in ENCODING_PROFILE["complete"]["encoding_names"][2]]
if len(filered_encodings_seq) > 0:
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

rule utils_collect_encodings:
    input:
         sequence_based_encodings_in=\
             expand(rules.meta_workflow_sequence_based_encodings.output,
                    normalized_dataset=ENCODING_PROFILE["all"]["datasets"]) \
                 if len(filered_encodings_seq) > 0 else [],
         structure_based_encodings_in=\
             expand(rules.meta_workflow_structure_based_encodings.output,
                    normalized_dataset=ENCODING_PROFILE["all"]["datasets"]) \
                 if len(filered_encodings_all) > 0 else [] + \
             expand(rules.meta_workflow_structure_based_encodings_windowed.output,
                    normalized_dataset=ENCODING_PROFILE["complete"]["datasets"]) \
                 if len(filered_encodings_complete) > 0 else []
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
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule predict:
