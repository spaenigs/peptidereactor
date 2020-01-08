import pandas as pd

DRY_RUN = config["dry_run"]
REDUCED = config["reduced"]
CORES = 8
DATASETS = ["neuropeptides"] #, "amps"]

def get_aaindex():
    df = pd.read_csv("apps/iFeature/data/AAindex.txt", sep="\t", index_col=0)
    df.columns = df.columns[1:].tolist() + ["NaN"]
    df = df.iloc[:, :-1]
    return df.index.to_list()

rule all:
    input:
         expand("data/{nds}/csv/aaindex/aaindex_{aaindex}.csv",
                nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
                aaindex=get_aaindex()),
         # expand("data/{nds}/csv/apaac/apaac_lambda_{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        lambda_val=list(range(1, 31))),
         # expand("data/{nds}/csv/cksaagp/cksaagp_gap_{gap_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        gap_val=range(1, 31)),
         # expand("data/{nds}/csv/cksaap/cksaap_gap_{gap_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        gap_val=range(1, 31)),
         # expand("data/{nds}/csv/eaac/eaac_window_{window_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        window_val=range(1, 31)),
         # expand("data/{nds}/csv/egaac/egaac_window_{window_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        window_val=range(1, 31)),
         # expand("data/{nds}/csv/geary/geary_nlag_{nlag_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        nlag_val=range(1, 31)),
         # expand("data/{nds}/csv/ksctriad/ksctriad_gap_{gap_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        gap_val=range(1, 31)),
         # expand("data/{nds}/csv/moran/moran_nlag_{nlag_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        nlag_val=range(1, 31)),
         # expand("data/{nds}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        nlag_val=range(1, 31)),
         # expand("data/{nds}/csv/paac/paac_lambda_{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        lambda_val=range(1, 31)),
         # expand("data/{nds}/csv/qsorder/qsorder_nlag_{nlag_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        nlag_val=range(1, 31)),
         # expand("data/{nds}/csv/socnumber/socnumber_nlag_{nlag_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        nlag_val=range(1, 31)),
         #
         # expand("data/{nds}/csv/psekraac_type1/psekraac_type1_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type2/psekraac_type2_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type3A/psekraac_type3A_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type3B/psekraac_type3B_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type4/psekraac_type4_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type5/psekraac_type5_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type6A/psekraac_type6A_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type6B/psekraac_type6B_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type6C/psekraac_type6C_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type7/psekraac_type7_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type8/psekraac_type8_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type9/psekraac_type9_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type10/psekraac_type10_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type11/psekraac_type11_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type12/psekraac_type12_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type13/psekraac_type13_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type14/psekraac_type14_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type15/psekraac_type15_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/psekraac_type16/psekraac_type16_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{nds}/csv/fft/fft_{aaindex}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        aaindex=get_aaindex()),
         # expand("data/{nds}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),
         #
         # expand("data/{nds}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        distance=[0,3,6,9,12]),
         # expand("data/{nds}/csv/delaunay/delaunay_{algorithm}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        algorithm=["average_distance", "total_distance", "cartesian_product",
         #                   "number_instances", "frequency_instances"]),
         # expand("data/{nds}/csv/waac/waac_{aaindex}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        aaindex=get_aaindex()),
         # expand("data/{nds}/csv/flgc/flgc_{aaindex}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        aaindex=get_aaindex()),
         # expand("data/{nds}/csv/fldpc/fldpc_{aaindex}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        aaindex=get_aaindex()),
         #
         # expand("data/{nds}/csv/ngram_a2/ngram_a2_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_a2/ngram_a2_lsv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_a2/ngram_a2_sv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{nds}/csv/ngram_a3/ngram_a3_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_a3/ngram_a3_lsv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_a3/ngram_a3_sv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{nds}/csv/ngram_e2/ngram_e2_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_e2/ngram_e2_lsv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_e2/ngram_e2_sv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{nds}/csv/ngram_e3/ngram_e3_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_e3/ngram_e3_lsv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_e3/ngram_e3_sv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{nds}/csv/ngram_s2/ngram_s2_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_s2/ngram_s2_lsv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_s2/ngram_s2_sv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{nds}/csv/ngram_s3/ngram_s3_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_s3/ngram_s3_lsv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{nds}/csv/ngram_s3/ngram_s3_sv_{dim}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{nds}/csv/distance_frequency/distance_frequency_dn_{nterminal}_dc_{cterminal}.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]),
         #        nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),
         #
         # expand("data/{nds}/csv/disorder.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/disorderb.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/disorderc.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/aac.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/binary.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/blosum62.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/ctdc.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/ctdd.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/ctdt.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/ctriad.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/dde.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/dpc.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/gaac.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/gdpc.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/gtpc.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/tpc.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/zscale.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/pssm.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/sseb.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/ssec.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/ta.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/asa.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/blomap.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/distance_distribution.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])),
         # expand("data/{nds}/csv/qsar.csv",
         #        nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"]))


rule create_datasets:
    input:
         "data/{ds}/seqs.fasta",
         "data/{ds}/classes.txt"
    output:
         "data/{ds}/sequence_length_distribution.svg",
         "data/{ds}_ds1/seqs.fasta",
         "data/{ds}_ds1/classes.txt",
         "data/{ds}_ds2/seqs.fasta",
         "data/{ds}_ds2/classes.txt"
    run:
         cmd = f"snakemake -s 01_create_datasets.smk --config dataset={{wildcards.ds}} cores={CORES}"
         if DRY_RUN:
            shell(cmd + " -n --quiet")
            for o in list(output):
                shell(f"touch {o}")
         else:
            shell(cmd)

rule reduce_datasets:
    input:
         "data/{nds}/seqs.fasta",
         "data/{nds}/classes.txt"
    output:
         "data/{nds}/seqs_reduced.fasta",
         "data/{nds}/classes_reduced.txt"
    run:
         from modlamp.core import read_fasta, save_fasta

         if REDUCED and not DRY_RUN:
             x = 5

             seqs, names = read_fasta(str(input[0]))
             save_fasta(str(output[0]), seqs[:x], names[:x])

             with open(str(input[1])) as f_in, open(str(output[1])) as f_out:
                 classes = list(map(lambda l: l.rstrip(), f_in.readlines()))
                 for c in classes[:x]:
                     f_out.write(f"{c}\n")
                     f_out.flush()

         else:
             shell("cp {input[0]} {output[0]}")
             shell("cp {input[1]} {output[1]}")

rule create_dirs:
    output:
         directory("data/{nds}/pdb/")

rule create_profiles:
    input:
         "data/{nds}/seqs_reduced.fasta",
         "data/{nds}/classes_reduced.txt",
         "data/{nds}/pdb/"
    output:
         "data/{nds}/seqs_msa.fasta",
         "data/{nds}/annotated_seqs.fasta",
         "data/{nds}/annotated_seqs_msa.fasta",
         "data/{nds}/annotated_classes.txt",
         "data/{nds}/annotated_pdbs_seqs.fasta",
         "data/{nds}/annotated_pdbs_classes.txt",
         directory("data/{nds}/profile/")
    run:
         cmd = f"snakemake -s 02_create_profiles.smk --config dataset={{wildcards.nds}} cores={CORES}"
         if DRY_RUN:
            shell(cmd + " -n --quiet")
            for o in list(output):
                shell(f"touch {o}")
         else:
            shell(cmd)

rule compute_params:
    input:
         "data/{nds}/annotated_seqs.fasta"
    output:
         expand("data/{{nds}}/misc/{encoding}.yaml",
                           encoding=["ksctriad", "moran", "nmbroto", "geary",
                                     "qsorder", "socnumber", "eaac", "cksaagp",
                                     "cksaap", "apaac", "paac"])
    run:
         cmd = f"snakemake -s 03_compute_params.smk --config dataset={{wildcards.nds}} cores={CORES}"
         if DRY_RUN:
            shell(cmd + " -n --quiet")
            for o in list(output):
                shell(f"touch {o}")
         else:
            shell(cmd)

rule create_encodings:
    input:
         "data/{nds}/seqs_msa.fasta",
         "data/{nds}/annotated_seqs.fasta",
         "data/{nds}/annotated_seqs_msa.fasta",
         "data/{nds}/annotated_classes.txt",
         "data/{nds}/annotated_pdbs_seqs.fasta",
         "data/{nds}/annotated_pdbs_classes.txt",
         "data/{nds}/profile/",
         expand("data/{{nds}}/misc/{encoding}.yaml",
                           encoding=["ksctriad", "moran", "nmbroto", "geary",
                                     "qsorder", "socnumber", "eaac", "cksaagp",
                                     "cksaap", "apaac", "paac"])
    output:
         expand("data/{{nds}}/csv/aaindex/aaindex_{aaindex}.csv",
                aaindex=get_aaindex()),
         # expand("data/{{nds}}/csv/apaac/apaac_lambda_{lambda_val}.csv",
         #        lambda_val=list(range(1, 31))),
         # expand("data/{{nds}}/csv/cksaagp/cksaagp_gap_{gap_val}.csv",
         #        gap_val=range(1, 31)),
         # expand("data/{{nds}}/csv/cksaap/cksaap_gap_{gap_val}.csv",
         #        gap_val=range(1, 31)),
         # expand("data/{{nds}}/csv/eaac/eaac_window_{window_val}.csv",
         #        window_val=range(1, 31)),
         # expand("data/{{nds}}/csv/egaac/egaac_window_{window_val}.csv",
         #        window_val=range(1, 31)),
         # expand("data/{{nds}}/csv/geary/geary_nlag_{nlag_val}.csv",
         #        nlag_val=range(1, 31)),
         # expand("data/{{nds}}/csv/ksctriad/ksctriad_gap_{gap_val}.csv",
         #        gap_val=range(1, 31)),
         # expand("data/{{nds}}/csv/moran/moran_nlag_{nlag_val}.csv",
         #        nlag_val=range(1, 31)),
         # expand("data/{{nds}}/csv/nmbroto/nmbroto_nlag_{nlag_val}.csv",
         #        nlag_val=range(1, 31)),
         # expand("data/{{nds}}/csv/paac/paac_lambda_{lambda_val}.csv",
         #        lambda_val=range(1, 31)),
         # expand("data/{{nds}}/csv/qsorder/qsorder_nlag_{nlag_val}.csv",
         #        nlag_val=range(1, 31)),
         # expand("data/{{nds}}/csv/socnumber/socnumber_nlag_{nlag_val}.csv",
         #        nlag_val=range(1, 31)),
         #
         # expand("data/{{nds}}/csv/psekraac_type1/psekraac_type1_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type2/psekraac_type2_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[2, 3, 4, 5, 6, 8, 15, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type3A/psekraac_type3A_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type3B/psekraac_type3B_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type4/psekraac_type4_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5, 8, 9, 11, 13, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type5/psekraac_type5_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[3, 4, 8, 10, 15, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type6A/psekraac_type6A_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 5, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type6B/psekraac_type6B_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type6C/psekraac_type6C_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[5],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type7/psekraac_type7_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type8/psekraac_type8_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type9/psekraac_type9_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type10/psekraac_type10_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type11/psekraac_type11_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type12/psekraac_type12_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,19)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type13/psekraac_type13_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=[4, 12, 17, 20],
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type14/psekraac_type14_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=list(range(2,21)),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type15/psekraac_type15_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/psekraac_type16/psekraac_type16_"
         #        "subtype-{sub_val}_raactype-{raac_val}_ktuple-{ktuple_val}_lambda-{lambda_val}.csv",
         #        sub_val=["g-gap", "lambda-correlation"], raac_val=(list(range(2,17)) + [20]),
         #        ktuple_val=list(range(1,4)), lambda_val=list(range(1,4))),
         # expand("data/{{nds}}/csv/fft/fft_{aaindex}.csv", aaindex=get_aaindex()),
         # expand("data/{{nds}}/csv/cgr/cgr_res_{resolution}_sf_{sfactor}.csv",
         #        resolution=[10, 20, 100, 200], sfactor=[0.5, 0.8632713]),
         #
         # expand("data/{{nds}}/csv/electrostatic_hull/electrostatic_hull_{distance}.csv",
         #        distance=[0,3,6,9,12]),
         # expand("data/{{nds}}/csv/delaunay/delaunay_{algorithm}.csv",
         #        algorithm=["average_distance", "total_distance", "cartesian_product",
         #                   "number_instances", "frequency_instances"]),
         # expand("data/{{nds}}/csv/waac/waac_{aaindex}.csv",
         #        aaindex=get_aaindex()),
         # expand("data/{{nds}}/csv/flgc/flgc_{aaindex}.csv",
         #        aaindex=get_aaindex()),
         # expand("data/{{nds}}/csv/fldpc/fldpc_{aaindex}.csv",
         #        aaindex=get_aaindex()),
         #
         # expand("data/{{nds}}/csv/ngram_a2/ngram_a2_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_a2/ngram_a2_lsv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_a2/ngram_a2_sv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{{nds}}/csv/ngram_a3/ngram_a3_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_a3/ngram_a3_lsv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_a3/ngram_a3_sv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{{nds}}/csv/ngram_e2/ngram_e2_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_e2/ngram_e2_lsv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_e2/ngram_e2_sv_{dim}.csv",
         #       dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{{nds}}/csv/ngram_e3/ngram_e3_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_e3/ngram_e3_lsv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_e3/ngram_e3_sv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{{nds}}/csv/ngram_s2/ngram_s2_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_s2/ngram_s2_lsv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_s2/ngram_s2_sv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{{nds}}/csv/ngram_s3/ngram_s3_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_s3/ngram_s3_lsv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         # expand("data/{{nds}}/csv/ngram_s3/ngram_s3_sv_{dim}.csv",
         #        dim=[1, 5, 20, 50, 100, 200, 300]),
         #
         # expand("data/{{nds}}/csv/distance_frequency/distance_frequency_dn_{nterminal}_dc_{cterminal}.csv",
         #        nterminal=[5, 10, 20, 50, 100], cterminal=[5, 10, 20, 50, 100]),
         #
         # "data/{nds}/csv/disorder.csv",
         # "data/{nds}/csv/disorderb.csv",
         # "data/{nds}/csv/disorderc.csv",
         # "data/{nds}/csv/aac.csv",
         # "data/{nds}/csv/binary.csv",
         # "data/{nds}/csv/blosum62.csv",
         # "data/{nds}/csv/ctdc.csv",
         # "data/{nds}/csv/ctdd.csv",
         # "data/{nds}/csv/ctdt.csv",
         # "data/{nds}/csv/ctriad.csv",
         # "data/{nds}/csv/dde.csv",
         # "data/{nds}/csv/dpc.csv",
         # "data/{nds}/csv/gaac.csv",
         # "data/{nds}/csv/gdpc.csv",
         # "data/{nds}/csv/gtpc.csv",
         # "data/{nds}/csv/tpc.csv",
         # "data/{nds}/csv/zscale.csv",
         # "data/{nds}/csv/pssm.csv",
         # "data/{nds}/csv/sseb.csv",
         # "data/{nds}/csv/ssec.csv",
         # "data/{nds}/csv/ta.csv",
         # "data/{nds}/csv/asa.csv",
         # "data/{nds}/csv/blomap.csv",
         # "data/{nds}/csv/distance_distribution.csv",
         # "data/{nds}/csv/qsar.csv"
    run:
         cmd = f"snakemake -s 04_create_encodings.smk --config dataset={{wildcards.nds}} cores={CORES}"
         if DRY_RUN:
            shell(cmd + " -n --quiet")
            for o in list(output):
                shell(f"touch {o}")
         else:
            shell(cmd)

if DRY_RUN:
    onsuccess:
        for d in expand("data/{nds}/", nds=expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])):
            shell(f"rm -rf {d}")