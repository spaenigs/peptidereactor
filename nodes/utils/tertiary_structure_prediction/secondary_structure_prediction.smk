from modlamp.core import read_fasta, save_fasta
import joblib as jl

import os
import textwrap

TOKEN = config["token"]
TARGET_PDBS = config["pdbs_out"]
TARGET_DIR = config["target_dir"]
PROFILE_DIR = config["profile_dir"]

PATH_TO_TMPS = "peptidereactor/RaptorX/tmp/"

rule get_psipred_profile:
    input:
         PROFILE_DIR + "{seq_name}.tgt"
    output:
         PROFILE_DIR + "{seq_name}.ss2",
         PROFILE_DIR + "{seq_name}.horiz"
    run:
         name = wildcards.seq_name
         shell(f"cp {PATH_TO_TMPS}{name}.ss2 {output[0]}")
         shell(f"cp {PATH_TO_TMPS}{name}.horiz {output[1]}")

rule get_vsl2_profile:
    input:
         PATH_TO_TMPS + f"{name}.diso"
    output:
         PROFILE_DIR + "{seq_name}.dis",
    run:
         name = wildcards.seq_name
         header = textwrap.dedent(
         """\
             Dummy file to mimic output from VSL2 Predictor of Intrinsically Disordered Regions
             
             Prediction Scores:
             ========================================
             NO.     RES.    PREDICTION      DISORDER
             ----------------------------------------
         """)
         with open(input[0]) as f, open(output[0]) as f2:
             for l in f.readlines()[5:]:
                 idx, aa, dis, prob, _  = l.rstrip().split("\n")
                 header += f"{idx}\t{aa}\t{100-prob}\t{'D' if dis == '*' else '.'}\n"
             footer = "========================================\n"
             f2.write(header + footer)
             f2.flush()

rule get_spx_profile:
    input:
         PATH_TO_TMPS + "tmp/{seq_name}.psp"
    output:
         PROFILE_DIR + "{seq_name}.spXout",
         temp(f"data/temp/{TOKEN}/protein-list-file_{{seq_name}}.txt")
    shell:
         f"""
         export spineXcodir=peptidereactor/spineXpublic;
         echo '{{wildcards.seq_name}}' > '{{output[1]}}';
         if [ -s {{input[0]}} ]
         then
             $spineXcodir/spX.pl '{{output[1]}}' {PROFILE_DIR} {PROFILE_DIR}
         else
             touch '{{output[0]}}';
             touch '{{output[1]}}';
         fi
         """

# TODO cut profiles for short sequences

rule remove_non_pssm_hits:
    input:
         config["fasta_in"],
         config["fasta_msa_in"],
         config["classes_in"],
         lambda wildcards: \
            expand(PROFILE_DIR + "{seq_name}.{ftype}",
                   ftype=["ss2", "horiz", "dis", "spXout"],
                   seq_name=read_fasta(config["fasta_in"])[1])
    output:
         temp(f"data/temp/{TOKEN}/filtered.joblib"),
         temp(f"data/temp/{TOKEN}/filtered_msa.joblib")
    run:
         def _filter(input_data_):
             filtered_seq_names = \
                 [fi.replace(PROFILE_DIR, "").replace(".ss2", "")
                  for fi in filter(lambda path: ".ss2" in path and os.stat(path).st_size > 0, list(input[3:]))]
             if len(filtered_seq_names) == 0:
                 return [], []
             else:
                 res_seqs, res_classes = zip(*filter(lambda tup: tup[0][0] in filtered_seq_names, zip(*input_data_)))
                 return res_seqs, res_classes

         seqs, names = read_fasta(str(input[0]))
         seqs_msa, names_msa = read_fasta(str(input[1]))
         with open(str(input[2])) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         in_da = [[[n, s] for n, s in zip(names, seqs)], classes]
         in_da_msa = [[[n, s] for n, s in zip(names_msa, seqs_msa)], classes]

         jl.dump(value=_filter(in_da), filename=str(output[0]))
         jl.dump(value=_filter(in_da_msa), filename=str(output[1]))

rule annotate_sequence_names:
    input:
         f"data/temp/{TOKEN}/filtered.joblib",
         f"data/temp/{TOKEN}/filtered_msa.joblib"
    output:
         config["fasta_sec_out"],
         config["fasta_sec_msa_out"],
         config["classes_sec"]
    run:
         def add_names(input_data_):
             res_seqs, res_classes = [], []
             for tup in zip(*input_data_):
                 seq_tup = tup[0]
                 seq_tup[0] = f"{seq_tup[0]}"
                 res_seqs.append(seq_tup)
                 res_classes.append(tup[1])
             return res_seqs, res_classes
         input_data, input_data_msa = jl.load(input[0], jl.load(input[1]))

         out = add_names(input_data)
         out_msa = add_names(input_data_msa)

         if len(out[0]) == 0:
             shell("touch {output[0]} {output[1]} {output[2]}")

         else:
             names, seqs = zip(*out[0])
             classes = out[1]
             save_fasta(output[0], seqs, names)

             names_msa, seqs_msa = zip(*out_msa[0])
             save_fasta(output[1], seqs_msa, names_msa)

             with open(output[2], mode="a") as f:
                 for c in out[1]:
                     f.write(f"{str(c)}\n")
                     f.flush()
