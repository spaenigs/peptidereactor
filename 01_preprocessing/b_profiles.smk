from sklearn.externals import joblib as jl
from Bio import SeqIO


localrules: split_input_data


rule split_input_data:
    input:
         ancient("00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed.joblib")
    output:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_{seq_name}_normal_distributed.joblib"
    group:
        "pssm"
    run:
        seq_tuple = list(filter(lambda t: wildcards.seq_name == t[0][0], zip(*jl.load(str(input)))))[0]
        jl.dump(value=([seq_tuple[0]], seq_tuple[1]), filename=str(output))


rule generate_pssm_profile:
    input:
        ancient("00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_{seq_name}_normal_distributed.joblib")
    output:
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.pssm",
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.asn.pssm",
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.mat"
    priority:
        1
    group:
        "pssm"
    run:
        import encoder.ifeature.pssm.utils as pssm_utils
        input_data = jl.load(str(input))
        input_data[0][0][0] = f"{wildcards.dataset}_{wildcards.part}_{wildcards.seq_name}"
        for o in output[:3]:
            with open(str(o), mode="w") as f:
                f.write("")
        pssm_utils.PSSMUtils.generate_profile(
            input_data,
            f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile",
            cores=1,
            db=config["uniref_db"])


rule generate_psipred_profile:
    input:
        ancient("00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_{seq_name}_normal_distributed.joblib"),
        ancient("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.asn.pssm")
    output:
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.ss2",
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.horiz",
        temp("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.ss"),
        temp("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.mtx")
    group:
        "pssm"
    shell:
        """
        export datadir={config[conda_env]}/share/psipred_4.01/data;
        if [ -s {input[1]} ]
        then
            chkparse '{input[1]}' > '{output[3]}';
            psipred '{output[3]}' $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 \
                > '{output[2]}';    
            psipass2 $datadir/weights_p2.dat 1 1.0 1.0 '{output[0]}' '{output[2]}' \
                > '{output[1]}';
        else
            touch '{output[0]}'
            touch '{output[1]}'
            touch '{output[2]}'
            touch '{output[3]}'
        fi
        """


rule generate_vsl2_profile:
    input:
        ancient("00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_{seq_name}_normal_distributed.joblib"),
        ancient("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.asn.pssm"),
        ancient("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.ss2")
    output:
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.dis",
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.flat"
    group:
        "pssm"
    run:
        input_data = jl.load(str(input[0]))
        seq_tup, class_ = input_data[0][0], input_data[1]
        fasta_name, fasta_seq = seq_tup[0], seq_tup[1]
        shell(f"echo {fasta_seq} > '{str(output[1])}'")
        shell("""
            if [ -s {input[1]} ] 
            then
                java -Duser.country=US -Duser.language=en -jar {config[cwd]}/{config[programs][vsl2]}/VSL2.jar -p:'{input[1]}' -s:'{output[1]}' -i:'{input[2]}' \
                    > '{output[0]}'
            else
                touch '{output[0]}'
                touch '{output[1]}'
            fi
        """)


rule generate_spx_profile:
    input:
        ancient("00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_{seq_name}_normal_distributed.joblib"),
        ancient("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.mat")
    output:
        "00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.spXout",
        temp("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_protein-list-file_{seq_name}.txt")
    group:
        "pssm"
    shell:
        """
        export spineXcodir={config[cwd]}/{config[programs][spinex]};
        echo '{wildcards.dataset}_{wildcards.part}_{wildcards.seq_name}' > '{output[1]}';
        if [ -s {input[1]} ]
        then
            {config[cwd]}/{config[programs][spinex]}/spX.pl '{output[1]}' 00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile 00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile
        else 
            touch '{output[0]}';
            touch '{output[1]}';
        fi
        """


rule generate_multiple_sequence_alignment:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed.joblib"
    output:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed_msa.joblib"
    run:
        from encoder.encoder import BaseEncoder as base_encoder
        input_data = jl.load(str(input))
        _, fastas_aligned = base_encoder.run_muscle(input_data)
        input_data = (fastas_aligned, input_data[1])
        jl.dump(value=input_data, filename=str(output))


def get_target_files(wildcards):
    fasta_file = checkpoints.save_as_fasta.get(dataset=config['dataset'], part=config['part']).output
    return expand("00_data/out/{dataset}/{dataset}_{part}/profile/{dataset}_{part}_{seq_name}.{ftype}",
                   dataset=config["dataset"], part=config["part"], ftype=["ss2", "horiz", "dis", "flat", "spXout", "mat", "pssm", "asn.pssm"],
                   seq_name=list(map(lambda r: str(r.name), SeqIO.parse(fasta_file[0], "fasta"))))

rule remove_non_pssm_hits:
    input:
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed.joblib",
         "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_normal_distributed_msa.joblib",
         get_target_files
    output:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib",
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered_msa.joblib"
    run:
        import os
        def _filter(input_data_):
            print(list(input[2:]))
            filtered_seq_names = \
                [fi.replace(f"00_data/out/{wildcards.dataset}/{wildcards.dataset}_{wildcards.part}/profile/{wildcards.dataset}_{wildcards.part}_", "")\
                     .replace(".mat", "")
                 for fi in filter(lambda path: ".mat" in path and os.stat(path).st_size > 0, list(input[2:]))]
            res_seqs, res_classes = zip(*filter(lambda tup: tup[0][0] in filtered_seq_names, zip(*input_data_)))
            return res_seqs, res_classes
        jl.dump(value=_filter(jl.load(str(input[0]))),
                filename=str(output[0]))
        jl.dump(value=_filter(jl.load(str(input[1]))),
                filename=str(output[1]))


rule annotate_sequence_names:
    input:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered.joblib",
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_pssms_filtered_msa.joblib"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated.joblib",
        "00_data/out/{dataset}/{dataset}_{part}/joblib/{dataset}_{part}_annotated_msa.joblib"
    run:
        def add_names(input_data_):
            res_seqs, res_classes = [], []
            for tup in zip(*input_data_):
                seq_tup = tup[0]
                print(f"original seq: {seq_tup[0]}")
                print(f"with wildcards: {wildcards.dataset}_{wildcards.part}_{seq_tup[0]}")
                seq_tup[0] = f"{wildcards.dataset}_{wildcards.part}_{seq_tup[0]}"
                res_seqs.append(seq_tup)
                res_classes.append(tup[1])
            return res_seqs, res_classes
        input_data, input_data_msa = jl.load(str(input[0])), jl.load(str(input[1]))
        jl.dump(value=add_names(input_data), filename=str(output[0]))
        jl.dump(value=add_names(input_data_msa), filename=str(output[1]))

