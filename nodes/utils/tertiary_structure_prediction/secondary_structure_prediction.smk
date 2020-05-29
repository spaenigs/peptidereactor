from modlamp.core import read_fasta, save_fasta

import pandas as pd
import joblib as jl

import os
import textwrap

TOKEN = config["token"]
TARGET_DIR = config["target_dir"]
PROFILE_DIR = config["profile_dir"]

PATH_TO_TMPS = "peptidereactor/RaptorX/tmp/"

rule get_vsl2_profile:
    input:
         PATH_TO_TMPS + "{seq_name}.diso"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.dis"
    run:
         if os.path.getsize(input[0]) != 0:
             name = wildcards.seq_name
             header = textwrap.dedent(
             """\
                 Dummy file to mimic output from VSL2 Predictor of Intrinsically Disordered Regions
                 
                 Prediction Scores:
                 ========================================
                 NO.     RES.    PREDICTION      DISORDER
                 ----------------------------------------
             """)
             with open(input[0]) as f, open(output[0], "w") as f2:
                 for l in f.readlines()[5:]:
                     idx, aa, dis, prob, _  = \
                         [i for i in l.lstrip().rstrip().split(" ") if i is not ""]
                     header += f"{idx}\t{aa}\t{100-float(prob)}\t{'D' if dis == '*' else '.'}\n"
                 footer = "========================================\n"
                 f2.write(header + footer)
                 f2.flush()
         else:
             shell("touch {output[0]}")

rule slice_or_dump_vsl2:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.dis",
         f"data/temp/{TOKEN}/{{seq_name}}.csv"
    output:
         PROFILE_DIR + "{seq_name}.dis"
    run:
         if os.path.getsize(input[0]) != 0:
             df = pd.read_csv(input[1])

             if df.empty:
                 shell("cp {input[0]} {output[0]}")

             else:
                 start, end = df["sstart"].values[0], df["send"].values[0]
                 with open(input[0]) as fi, open(output[0], "w") as fo:
                     lines = fi.readlines()
                     header = lines[:6]
                     footer = lines[-1:]
                     for l in lines[6:]:
                         if l.startswith("="):
                             continue
                         line_splitted = l.lstrip().rstrip().split("\t")
                         id = int(line_splitted[0])
                         if start <= id <= end:
                             header += [l]
                     for h in header + footer:
                         fo.write(h)
                         fo.flush()
         else:
             shell("touch {output[0]}")

rule get_spx_profile:
    input:
         PATH_TO_TMPS + "{seq_name}.psp",
         "peptidereactor/spineXpublic/spineXpublic_moved.txt"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.spXout",
         temp(f"data/temp/{TOKEN}/protein-list-file_{{seq_name}}.txt")
    run:
         if os.path.getsize(input[0]) != 0:
             shell(
             f"""
             cp {{input[0]}} data/temp/{TOKEN}/{{wildcards.seq_name}}.mat;
             export spineXcodir=peptidereactor/spineXpublic;
             echo '{{wildcards.seq_name}}' > '{{output[1]}}';
             if [ -s {{input[0]}} ]
             then
                 $spineXcodir/spX.pl '{{output[1]}}' data/temp/{TOKEN}/ data/temp/{TOKEN}/
             else
                 touch '{{output[0]}}';
                 touch '{{output[1]}}';
             fi
             """)
         else:
             shell("touch {output[0]} {output[1]}")

rule slice_or_dump_spx:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.spXout",
         f"data/temp/{TOKEN}/{{seq_name}}.csv"
    output:
         PROFILE_DIR + "{seq_name}.spXout"
    run:
         if os.path.getsize(input[0]) != 0:

             df = pd.read_csv(input[1])

             if df.empty:
                 shell("cp {input[0]} {output[0]}")

             else:
                 start, end = df["sstart"].values[0], df["send"].values[0]
                 with open(input[0]) as fi, open(output[0], "w") as fo:
                     lines = fi.readlines()
                     header = lines[:1]
                     for l in lines[1:]:
                         line_splitted = l.lstrip().rstrip().split(" ")
                         id = int(line_splitted[0])
                         if start <= id <= end:
                             header += [l]
                     for h in header:
                         fo.write(h)
                         fo.flush()
         else:
             shell("touch {output[0]}")

rule remove_non_pssm_hits:
    input:
         config["fasta_in"],
         config["classes_in"],
         lambda wildcards: \
            expand(PROFILE_DIR + "{seq_name}.{ftype}", ftype=["dis", "spXout"],
                   seq_name=read_fasta(config["fasta_in"])[1])
    output:
         temp(f"data/temp/{TOKEN}/filtered.joblib"),
    run:
         def _filter(input_data_):
             filtered_seq_names = \
                 [fi.replace(PROFILE_DIR, "").replace(".spXout", "")
                  for fi in filter(lambda path: ".spXout" in path and os.stat(path).st_size > 0, input[2:])]
             if len(filtered_seq_names) == 0:
                 return [], []
             else:
                 res_seqs, res_classes = zip(*filter(lambda tup: tup[0][0] in filtered_seq_names, zip(*input_data_)))
                 return res_seqs, res_classes

         seqs, names = read_fasta(input[0])
         with open(input[1]) as f:
             classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         in_da = [[[n, s] for n, s in zip(names, seqs)], classes]
         jl.dump(value=_filter(in_da), filename=output[0])

rule annotate_sequence_names:
    input:
         f"data/temp/{TOKEN}/filtered.joblib",
    output:
         config["fasta_sec_out"],
         config["classes_sec_out"]
    run:
         def replace_seqs(input_data_):
             res_seqs, res_classes = [], []
             for tup in zip(*input_data_):
                 seq_tup = tup[0]
                 df = pd.read_csv(f"data/temp/{TOKEN}/{seq_tup[0]}.csv", engine="c")
                 if not df.empty:
                    seq_tup[1] = df["sseq"][0]
                 res_seqs.append(seq_tup)
                 res_classes.append(tup[1])
             return res_seqs, res_classes

         input_data = jl.load(input[0])
         out = replace_seqs(input_data)

         if len(out[0]) == 0:
             shell("touch {output[0]} {output[1]}")

         else:
             names, seqs = zip(*out[0])
             classes = out[1]
             save_fasta(output[0], seqs, names)

             with open(output[1], mode="a") as f:
                 for c in out[1]:
                     f.write(f"{str(c)}\n")
                     f.flush()
