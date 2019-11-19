TOKEN = config["token"]

rule init:
    input:
         config["download_link_in"]
    output:
         "apps/RaptorX/setup.pl"
    priority:
         50
    shell:
         f"""
         export link=$(head -n 1 {{input}}); 
         wget $link -P data/temp/{TOKEN}/;
         wget $(cat data/temp/{TOKEN}/index.html | grep CNFsearch | \
            perl -ne 'print "$1\n" if /(http.*?)\W>/') -O data/temp/{TOKEN}/CNFsearch.zip \
            -q --show-progress --progress=bar:force:noscroll;
         unzip -q data/temp/{TOKEN}/CNFsearch.zip -d data/temp/{TOKEN}/; 
         mv data/temp/{TOKEN}/CNFsearch1.66_complete/* apps/RaptorX/;   
         export OLDWD=$PWD; cd apps/RaptorX/; ./setup.pl; cd $OLDWD;
         """

rule set_modeller_license_key:
    input:
         config["license_key_in"]
    output:
         f"data/temp/{TOKEN}/config.py"
    shell:
         """
         export key=$(head -n 1 {input});
         export config_file=$(mod9.23 -v 2> >(grep config.py));
         perl -i -pe "s/XXXX/$key/g" $config_file;
         export stderr_out=$(mod9.23 -v 2> >(grep 'Cannot open file -v'))
         if [ $stderr_out = "" ]; then 
            echo "Check your MODELLER license key.";
         else
            cp $config_file {output}
         fi
         """

rule download_databases:
    input:
         ancient(config["download_link_in"]),
         ancient("apps/RaptorX/setup.pl")
    output:
         "apps/RaptorX/databases/{target_folder}/{database}{type}.tar.gz"
    shell:
         f"""
         export link=$(head -n 1 {{input[0]}}); 
         wget $link -P data/temp/{TOKEN}/;
         wget $(cat data/temp/{TOKEN}/index.html | grep nr{{wildcards.size}}.tar.gz | \
            perl -ne 'print "$1\n" if /(http.*?)\W>/') -O {{output}} \
            -q --show-progress --progress=bar:force:noscroll; 
         """

rule collect:
    input:
         expand("apps/RaptorX/databases/{folder}/{database}{type}.tar.gz",
                folder="")

rule download_nr:
    input:
         ancient(config["download_link_in"]),
         ancient("apps/RaptorX/setup.pl")
    output:
         "apps/RaptorX/databases/NR_new/nr{size}.tar.gz"
    priority:
         50
    shell:
         f"""
         export link=$(head -n 1 {{input[0]}}); 
         wget $link -P data/temp/{TOKEN}/;
         wget $(cat data/temp/{TOKEN}/index.html | grep nr{{wildcards.size}}.tar.gz | \
            perl -ne 'print "$1\n" if /(http.*?)\W>/') -O {{output}} \
            -q --show-progress --progress=bar:force:noscroll; 
         """

rule unzip_nr:
    input:
         "apps/RaptorX/databases/NR_new/nr{size}.tar.gz"
    output:
         "apps/RaptorX/databases/NR_new/nr{size}.pal"
    priority:
         30
    shell:
         f"""
         tar -zxf {{input}} -C data/temp/{TOKEN}/;
         mv data/temp/{TOKEN}/nr{{wildcards.size}}* apps/RaptorX/databases/NR_new/;
         """

rule download_tpl:
    input:
         ancient(config["download_link_in"]),
         ancient("apps/RaptorX/setup.pl")
    output:
         "apps/RaptorX/databases/TPL_{type}/TPL_{type}.tar.gz"
    priority:
         50
    shell:
         f"""
         export link=$(head -n 1 {{input[0]}}); 
         wget $link -P data/temp/{TOKEN}/;
         wget $(cat data/temp/{TOKEN}/index.html | grep TPL_{{wildcards.type}}.*.tar.gz | \
            perl -ne 'print "$1\n" if /(http.*?)\W>/') -O {{output}} \
            -q --show-progress --progress=bar:force:noscroll;
         """

rule unzip_tpl:
    input:
         "apps/RaptorX/databases/TPL_{type}/TPL_{type}.tar.gz"
    output:
         "apps/RaptorX/databases/TPL_{type}/TPL_{type}.finished.txt"
    priority:
         30
    shell:
         f"""
         tar -zxf {{input}} -C data/temp/{TOKEN}/;
         mv data/temp/{TOKEN}/TPL_{{wildcards.type}}/* apps/RaptorX/databases/TPL_{{wildcards.type}}/
         touch {{output}};
         """

rule move_tpl_remain:
    input:
         "apps/RaptorX/databases/TPL_Remain/TPL_Remain.finished.txt"
    output:
         "apps/RaptorX/databases/TPL_BC100/TPL_BC100.finished.txt"
    shell:
         """
         for file in `ls -1 apps/RaptorX/databases/TPL_Remain/`; do 
            mv apps/RaptorX/databases/TPL_Remain/$file apps/RaptorX/databases/TPL_BC100/$file; 
         done
         rm -r apps/RaptorX/databases/TPL_Remain/;
         touch {output};
         """

rule download_tpl_list:
    input:
         ""
    output:
         ""
    run:
         ""

rule download_pdb:
    input:
         ancient(config["download_link_in"]),
         ancient("apps/RaptorX/setup.pl")
    output:
         "apps/RaptorX/databases/pdb_{type}/pdb_{type}.tar.gz"
    priority:
         50
    shell:
         f"""
         export link=$(head -n 1 {{input[0]}}); 
         wget $link -P data/temp/{TOKEN}/;
         wget $(cat data/temp/{TOKEN}/index.html | grep pdb_{{wildcards.type}}.*.tar.gz | \
            perl -ne 'print "$1\n" if /(http.*?)\W>/') -O {{output}} \
            -q --show-progress --progress=bar:force:noscroll;
         """

rule unzip_pdb:
    input:
         "apps/RaptorX/databases/pdb_{type}/pdb_{type}.tar.gz"
    output:
         "apps/RaptorX/databases/pdb_{type}/pdb_{type}.finished.txt"
    priority:
         30
    shell:
         f"""
         tar -zxf {{input}} -C data/temp/{TOKEN}/;
         mv data/temp/{TOKEN}/pdb_{{wildcards.type}}/* apps/RaptorX/databases/pdb_{{wildcards.type}}/
         touch {{output}};
         """

rule move_pdb_remain:
    input:
         "apps/RaptorX/databases/pdb_Remain/pdb_Remain.finished.txt"
    output:
         "apps/RaptorX/databases/pdb_BC100/pdb_BC100.finished.txt"
    shell:
         """
         for file in `ls -1 apps/RaptorX/databases/pdb_Remain/`; do 
            mv apps/RaptorX/databases/pdb_Remain/$file apps/RaptorX/databases/pdb_BC100/$file; 
         done
         rm -r apps/RaptorX/databases/pdb_Remain/;
         touch {output};
         """

# def generate_pdb_file_names(wildcards):
#     files = []
#     for ds in config["datasets"]["processed"]["names"]:
#         files += expand("snakemake/pdb/{dataset}-{seq_name}.pdb",
#                         dataset=ds,
#                         seq_name=[i.rstrip().replace(">", "") for i in open(f"snakemake/fasta/{ds}.fasta").readlines()[0::2]])
#     return files