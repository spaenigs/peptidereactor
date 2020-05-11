import os

from nodes.utils.tertiary_structure_search.scripts.utils \
    import get_seq_names

TOKEN = config["token"]
TARGET_PDBS = config["pdbs_out"]
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
         with open(str(input)) as f, open(str(output), mode="w") as out:
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
         TARGET_DIR + "{seq_name}.pdb"
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

# TODO cut profiles for short sequences

rule remove_non_hits:
    input:
         config["fasta_in"],
         config["classes_in"],
         expand(TARGET_DIR + "{seq_name}.pdb",
                seq_name=get_seq_names(config["fasta_in"]))
    output:
         config["fasta_out"],
         config["classes_out"]
    run:
         from modlamp.core import save_fasta, read_fasta

         seqs, names = read_fasta(str(input[0]))
         with open(str(input[1])) as f:
            classes = list(map(lambda l: int(l.rstrip()), f.readlines()))
         in_da = [[[n, s] for n, s in zip(names, seqs)], classes]

         seq_tuples = dict((name, tup) for name, tup in zip(names, zip(seqs, classes)))
         res_names, res_seqs, res_classes = [], [], []
         for file_path in list(input[2:]):
             if os.path.getsize(file_path) > 0:
                 seq_name = os.path.basename(file_path).replace(".pdb", "")
                 seq, class_ = seq_tuples[seq_name]
                 res_names.append(seq_name)
                 res_seqs.append(seq)
                 res_classes.append(class_)

         save_fasta(str(output[0]), sequences=res_seqs, names=res_names)

         with open(str(output[1]), mode="a") as f:
             for c in res_classes:
                 f.write(f"{c}\n")
