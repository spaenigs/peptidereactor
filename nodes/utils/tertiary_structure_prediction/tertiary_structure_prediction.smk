from Bio.PDB import Dice, PDBParser
from modlamp.core import save_fasta, read_fasta

import pandas as pd

import os

TOKEN = config["token"]
TARGET_DIR = config["target_dir"]
PROFILE_DIR = config["profile_dir"]

rule search_template_database:
    input:
         PROFILE_DIR + "{seq_name}.tgt",
         expand("peptidereactor/RaptorX/databases/{database}/{database}.moved.txt",
                database=["TPL_BC40", "TPL_Remain"])
    output:
         PROFILE_DIR + "{seq_name}.rank"
    shell:
         f"""
         touch {{output}};
         if [ -s {{input[0]}} ]; then
            export OLDWD=$PWD; cd peptidereactor/RaptorX/;
            ./CNFsearch -a 1 -q {{wildcards.seq_name}} -g $OLDWD/{PROFILE_DIR} -o $OLDWD/{{output}}
            cd - 1> /dev/null;
         else
            touch {{output}};
         fi
         """

rule best_template:
    input:
         PROFILE_DIR + "{seq_name}.rank"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}_best_template.txt")
    run:
         import re
         with open(input[0]) as f, open(output[0], mode="w") as out:
             lines_with_hit = list(filter(lambda line: re.search("^1\s", line) is not None,
                                          f.readlines()))
             if len(lines_with_hit) == 1:
                 line_with_hit = lines_with_hit[0]
                 id = line_with_hit.split()[1]
             else:
                 id = ""
             out.write(id)
             out.flush()

rule align_to_template:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}_best_template.txt"
    output:
         temp(f"data/temp/{TOKEN}/{{seq_name}}.raptorx.fasta"),
         temp(f"data/temp/{TOKEN}/{{seq_name}}.raptorx.cnfpred")
    shell:
         f"""
         if [ -s {{input}} ]; then
            export BEST_HIT=`head -n 1 {{input}}`;
            export OLDWD=$PWD; cd peptidereactor/RaptorX/; 
            ./CNFalign_lite -q {{wildcards.seq_name}} \
                            -t $BEST_HIT \
                            -l databases/TPL_BC100/ \
                            -g $OLDWD/{PROFILE_DIR} \
                            -d $OLDWD/data/temp/{TOKEN}/;
            cd - 1> /dev/null;
            mv data/temp/{TOKEN}/$BEST_HIT-{{wildcards.seq_name}}.fasta {{output[0]}};
            mv data/temp/{TOKEN}/$BEST_HIT-{{wildcards.seq_name}}.cnfpred {{output[1]}};
         else
            touch {{output[0]}}; touch {{output[1]}};
         fi         
         """

rule generate_structure:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.raptorx.fasta",
         f"data/temp/{TOKEN}/{{seq_name}}.raptorx.cnfpred",
         expand("peptidereactor/RaptorX/databases/{database}/{database}.moved.txt",
                database=["pdb_BC40", "pdb_Remain"]),
         "peptidereactor/RaptorX/modeller/config.py"
    output:
         f"data/temp/{TOKEN}/{{seq_name}}.pdbtmp",
    shell:
         """
         if [ -s {input[0]} ]; then
            export base_path=$(cat {input[4]} | \
                grep install_dir | \
                perl -ne 'print $1 if /.*? r\W([\-\.\d\w\/]*)/')
            export full_path=$base_path/modlib/modeller/
            cp {input[4]} $full_path;  
            export OLDWD=$PWD; cd peptidereactor/RaptorX/;
            ./build3Dmodel -i $OLDWD/{input[0]} \
                           -q {wildcards.seq_name} \
                           -d databases/pdb_BC100/ \
                           -m mod9.23;
            mv *{wildcards.seq_name}.*.pdb $OLDWD/{output};
            cd - 1> /dev/null;
         else 
            touch {output};
         fi
         """

rule slice_or_dump_ter:
    input:
         f"data/temp/{TOKEN}/{{seq_name}}.csv",
         f"data/temp/{TOKEN}/{{seq_name}}.pdbtmp"
    output:
         TARGET_DIR + "{seq_name}.pdb"
    run:
         df = pd.read_csv(input[0])

         if df.empty:
             shell("cp {input[1]} {output[0]}")

         elif os.path.getsize(input[1]) == 0:
             shell("cp {input[1]} {output[0]}")

         else:
             structure = PDBParser().get_structure(wildcards.seq_name, input[1])
             chain_id = [c.get_id() for c in structure.get_chains()][0]
             start, end = df["sstart"].values[0], df["send"].values[0]
             Dice.extract(structure, chain_id, start, end, output[0])

rule remove_non_hits:
    input:
         config["fasta_in"],
         config["classes_in"],
         expand(TARGET_DIR + "{seq_name}.pdb",
                seq_name=read_fasta(config["fasta_in"])[1])
    output:
         config["fasta_ter_out"],
         config["classes_ter_out"]
    run:
         seqs, names = read_fasta(input[0])
         with open(str(input[1])) as f:
            classes = list(map(lambda l: int(l.rstrip()), f.readlines()))

         in_da = [[[n, s] for n, s in zip(names, seqs)], classes]

         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         res_names, res_seqs, res_classes = [], [], []
         for file_path in input[2:]:
             if os.path.getsize(file_path) > 0:
                 seq_name = os.path.basename(file_path).replace(".pdb", "")
                 seq, class_ = seq_tuples[seq_name]
                 df = pd.read_csv(f"data/temp/{TOKEN}/{seq_name}.csv", engine="c")
                 if not df.empty:
                    seq = df["sseq"][0]
                 res_names.append(seq_name)
                 res_seqs.append(seq)
                 res_classes.append(class_)

         save_fasta(str(output[0]), sequences=res_seqs, names=res_names)

         with open(output[1], mode="a") as f:
             for c in res_classes:
                 f.write(f"{c}\n")
