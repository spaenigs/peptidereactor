from modlamp.core import read_fasta
import pandas as pd
import yaml

DATASETS = ["neuropeptides"]
CORES = 32

# checkpoint create_all_datasets:
#     input:
#          "data/{dataset}/seqs.fasta",
#          "data/{dataset}/classes.txt"
#     output:
#          expand("data/{dataset}_ds1/seqs.fasta", dataset=DATASETS),
#          expand("data/{dataset}_ds1/classes.txt", dataset=DATASETS),
#          expand("data/{dataset}_ds2/seqs.fasta", dataset=DATASETS),
#          expand("data/{dataset}_ds2/classes.txt", dataset=DATASETS)
#     shell:
#          f"./run_pipeline -s 01_create_dataset.smk --config dataset={{wildcards.dataset}} cores={CORES} --quiet"

rule all:
    input:
         expand("data/{dataset}_{id}/csv/disorder.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/disorderb.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/disorderc.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/aac.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/binary.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/blosum62.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/ctdc.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/ctdd.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/ctdt.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/ctriad.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/dde.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/dpc.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/gaac.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/gdpc.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/gtpc.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/tpc.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/zscale.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/pssm.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/sseb.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/ssec.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/ta.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/asa.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/blomap.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/distance_distribution.csv",
                dataset=DATASETS, id=["ds1", "ds2"]),
         expand("data/{dataset}_{id}/csv/qsar.csv",
                dataset=DATASETS, id=["ds1", "ds2"])

rule create_all_profiles:
    input:
         "data/{dataset}_{id}/seqs.fasta",
         "data/{dataset}_{id}/classes.txt"
    output:
         "data/{dataset}_{id}/seqs_msa.fasta",
         "data/{dataset}_{id}/annotated_seqs.fasta",
         "data/{dataset}_{id}/annotated_seqs_msa.fasta",
         "data/{dataset}_{id}/annotated_classes.txt",
         expand("data/{dataset}_{id}/profile/{{seq_name}}.{{ftype}}",
                seq_name=read_fasta("data/{dataset}_ds{id}/seqs.fasta")[1],
                ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"]),
         "data/{dataset}_ds{id}/annotated_pdbs_seqs.fasta",
         "data/{dataset}_ds{id}/annotated_pdbs_classes.txt",
         expand("data/{dataset}_ds{id}/pdb/{{seq_name}}.pdb",
                seq_name=read_fasta("data/{dataset}_ds{id}/seqs.fasta")[1])
    shell:
         f"./run_pipeline -s 02_create_profiles.smk --config dataset={{wildcards.dataset}} cores={CORES} --quiet"

rule compute_all_params:
    input:
         "data/{dataset}_ds{id}/annotated_seqs.fasta"
    output:
          expand("data/{dataset}_ds{id}/misc/{{encoding}}.yaml",
                 encoding=["ksctriad", "moran", "nmbroto", "geary",
                           "qsorder", "socnumber", "eaac", "cksaagp",
                           "cksaap", "apaac", "paac"]),
          expand("data/{dataset}_ds{id}/misc/ngram_{{type}}{{size}}.yaml",
                 type=["a","e","s"], size=[2,3])
    shell:
          f"./run_pipeline -s 03_compute_params.smk --config dataset={{wildcards.dataset}} cores={CORES} --quiet"

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

def get_max_vals(encoding, dataset):
    try:
        with open(f"data/{dataset}/misc/{encoding}.yaml") as f:
            return yaml.safe_load(f) + 1  # range is exclusive
    except FileNotFoundError:
        exit("""
        Please run node window_length beforehand (or set them manually).
        See, e.g., apps/iFeature/codes/KSCTriad.py for details.
        """)

def get_max_dim_size(ngram, dataset):
    ngram_type, size = list(ngram)
    try:
        with open(f"data/{dataset}/misc/ngram_{ngram_type}{size}.yaml") as f:
            return yaml.safe_load(f)  # range is exclusive
    except FileNotFoundError:
        exit("""
        Please run node dim_size beforehand (or set dimension manually): min(len(shape[0], shape[1]).
        """)

rule create_all_encodings:
    input:
         "data/{dataset}_ds{id}/seqs_msa.fasta",
         "data/{dataset}_ds{id}/annotated_seqs.fasta",
         "data/{dataset}_ds{id}/annotated_seqs_msa.fasta",
         "data/{dataset}_ds{id}/annotated_classes.txt",
         expand("data/{dataset}_ds{id}/profile/{{seq_name}}.{{ftype}}",
                seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1],
                ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"]),
         "data/{dataset}_ds{id}/annotated_pdbs_seqs.fasta",
         "data/{dataset}_ds{id}/annotated_pdbs_classes.txt",
         expand("data/{dataset}_ds{id}/pdb/{{seq_name}}.pdb",
                seq_name=read_fasta(f"data/{config['dataset']}/seqs.fasta")[1]),
         expand("data/{dataset}_ds{id}/misc/{{encoding}}.yaml",
                 encoding=["ksctriad", "moran", "nmbroto", "geary",
                           "qsorder", "socnumber", "eaac", "cksaagp",
                           "cksaap", "apaac", "paac"]),
         expand("data/{dataset}_ds{id}/misc/ngram_{{type}}{{size}}.yaml",
                type=["a","e","s"], size=[2,3])
    output:
         expand("data/{dataset}_ds{id}/csv/aaindex/aaindex_{{aaindex}}.csv",
               aaindex=get_aaindex()),
         expand("data/{dataset}_ds{id}/csv/apaac/apaac_lambda_{{lambda_val}}.csv",
                lambda_val=list(range(1, get_max_vals("apaac", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/cksaagp/cksaagp_gap_{{gap_val}}.csv",
                gap_val=list(range(1, get_max_vals("cksaagp", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/cksaap/cksaap_gap_{{gap_val}}.csv",
                gap_val=list(range(1, get_max_vals("cksaap", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/eaac/eaac_window_{{window_val}}.csv",
               window_val=list(range(1, get_max_vals("eaac", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/egaac/egaac_window_{{window_val}}.csv",
               window_val=list(range(1, 31))),
         expand("data/{dataset}_ds{id}/csv/geary/geary_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("geary", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/ksctriad/ksctriad_gap_{{gap_val}}.csv",
               gap_val=list(range(1, get_max_vals("ksctriad", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/moran/moran_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("moran", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/nmbroto/nmbroto_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("nmbroto", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/paac/paac_lambda_{{lambda_val}}.csv",
               lambda_val=list(range(1, get_max_vals("paac", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/qsorder/qsorder_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("qsorder", "{dataset}_ds{id}")))),
         expand("data/{dataset}_ds{id}/csv/socnumber/socnumber_nlag_{{nlag_val}}.csv",
               nlag_val=list(range(1, get_max_vals("socnumber", "{dataset}_ds{id}")))),

         expand("data/{dataset}_ds{id}/csv/psekraac_type1/psekraac_type1_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type2/psekraac_type2_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type3A/psekraac_type3A_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type3B/psekraac_type3B_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type4/psekraac_type4_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type5/psekraac_type5_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type6A/psekraac_type6A_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type6B/psekraac_type6B_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type6C/psekraac_type6C_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type7/psekraac_type7_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type8/psekraac_type8_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type9/psekraac_type9_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type10/psekraac_type10_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type11/psekraac_type11_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type12/psekraac_type12_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type13/psekraac_type13_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type14/psekraac_type14_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type15/psekraac_type15_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/psekraac_type16/psekraac_type16_"
               "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
               sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
               ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         expand("data/{dataset}_ds{id}/csv/fft/fft_{{aaindex}}.csv", aaindex=get_aaindex()),
         expand("data/{dataset}_ds{id}/csv/cgr/cgr_res_{{resolution}}_sf_{{sfactor}}.csv",
               resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),

         expand("data/{dataset}_ds{id}/csv/electrostatic_hull/electrostatic_hull_{{distance}}.csv",
               distance=[0,3,6,9,12]),
         expand("data/{dataset}_ds{id}/csv/delaunay/delaunay_{{algorithm}}.csv",
               algorithm=["average_distance", "total_distance", "cartesian_product",
                          "number_instances", "frequency_instances"]),
         expand("data/{dataset}_ds{id}/csv/waac/waac_{{aaindex}}.csv",
               aaindex=get_aaindex()),
         expand("data/{dataset}_ds{id}/csv/flgc/flgc_{{aaindex}}.csv",
               aaindex=get_aaindex()),
         expand("data/{dataset}_ds{id}/csv/fldpc/fldpc_{{aaindex}}.csv",
               aaindex=get_aaindex()),

         expand("data/{dataset}_ds{id}/csv/ngram_a2/ngram_a2_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a2", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_a2/ngram_a2_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a2", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_a2/ngram_a2_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a2", "{dataset}_ds{id}"))),

         expand("data/{dataset}_ds{id}/csv/ngram_a3/ngram_a3_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a3", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_a3/ngram_a3_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a3", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_a3/ngram_a3_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("a3", "{dataset}_ds{id}"))),

         expand("data/{dataset}_ds{id}/csv/ngram_e2/ngram_e2_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e2", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_e2/ngram_e2_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e2", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_e2/ngram_e2_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e2", "{dataset}_ds{id}"))),

         expand("data/{dataset}_ds{id}/csv/ngram_e3/ngram_e3_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e3", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_e3/ngram_e3_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e3", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_e3/ngram_e3_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("e3", "{dataset}_ds{id}"))),

         expand("data/{dataset}_ds{id}/csv/ngram_s2/ngram_s2_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s2", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_s2/ngram_s2_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s2", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_s2/ngram_s2_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s2", "{dataset}_ds{id}"))),

         expand("data/{dataset}_ds{id}/csv/ngram_s3/ngram_s3_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s3", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_s3/ngram_s3_lsv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s3", "{dataset}_ds{id}"))),
         expand("data/{dataset}_ds{id}/csv/ngram_s3/ngram_s3_sv_{{dim}}.csv",
               dim=range(1, get_max_dim_size("s3", "{dataset}_ds{id}"))),

         expand("data/{dataset}_ds{id}/csv/distance_frequency/distance_frequency_dn_{{nterminal}}_dc_{{cterminal}}.csv",
               nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),

         "data/{dataset}_ds{id}/csv/disorder.csv",
         "data/{dataset}_ds{id}/csv/disorderb.csv",
         "data/{dataset}_ds{id}/csv/disorderc.csv",
         "data/{dataset}_ds{id}/csv/aac.csv",
         "data/{dataset}_ds{id}/csv/binary.csv",
         "data/{dataset}_ds{id}/csv/blosum62.csv",
         "data/{dataset}_ds{id}/csv/ctdc.csv",
         "data/{dataset}_ds{id}/csv/ctdd.csv",
         "data/{dataset}_ds{id}/csv/ctdt.csv",
         "data/{dataset}_ds{id}/csv/ctriad.csv",
         "data/{dataset}_ds{id}/csv/dde.csv",
         "data/{dataset}_ds{id}/csv/dpc.csv",
         "data/{dataset}_ds{id}/csv/gaac.csv",
         "data/{dataset}_ds{id}/csv/gdpc.csv",
         "data/{dataset}_ds{id}/csv/gtpc.csv",
         "data/{dataset}_ds{id}/csv/tpc.csv",
         "data/{dataset}_ds{id}/csv/zscale.csv",
         "data/{dataset}_ds{id}/csv/pssm.csv",
         "data/{dataset}_ds{id}/csv/sseb.csv",
         "data/{dataset}_ds{id}/csv/ssec.csv",
         "data/{dataset}_ds{id}/csv/ta.csv",
         "data/{dataset}_ds{id}/csv/asa.csv",
         "data/{dataset}_ds{id}/csv/blomap.csv",
         "data/{dataset}_ds{id}/csv/distance_distribution.csv",
         "data/{dataset}_ds{id}/csv/qsar.csv"


