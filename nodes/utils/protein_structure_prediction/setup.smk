TOKEN = config["token"]

rule all:
    input:
        expand(f"data/temp/{TOKEN}/{{database}}_download_link.txt",
               database=["nr70", "nr90", "TPL_BC40", "TPL_Remain", "pdb_BC40", "pdb_Remain", "CNFsearch"])

rule download_index_html:
    input:
         config["download_link_in"]
    output:
         temp(f"data/temp/{TOKEN}/index.html")
    run:
         with open(str(input)) as file_in:
              link_to_index_html = list(file_in.readlines())[0]
              shell(f"wget {link_to_index_html} -P data/temp/{TOKEN}/ 2> /dev/null")

rule parse_download_links:
    input:
         f"data/temp/{TOKEN}/index.html"
    output:
         temp(f"data/temp/{TOKEN}/{{database}}_download_link.txt")
    run:
         import re
         with open(str(input)) as index_html:
             for line in index_html.readlines():
                 if re.match(f".*?({wildcards.database}).*", line) is not None:
                     download_link = re.match(".*href='(.*?)'", line).group(1)
         with open(str(output), mode="w") as file_out:
             file_out.write(download_link)
             file_out.flush()

rule init_raptorx:
    input:
         f"data/temp/{TOKEN}/CNFsearch_download_link.txt"
    output:
         "apps/RaptorX/setup.pl",
         temp(f"data/temp/{TOKEN}/CNFsearch.zip")
    priority:
         50
    shell:
         f"""
         wget $(head -n 1 {{input}}) \
            -O data/temp/{TOKEN}/CNFsearch.zip \
            -q --show-progress --progress=bar:force:noscroll;   
         unzip -q data/temp/{TOKEN}/CNFsearch.zip -d data/temp/{TOKEN}/; 
         mv data/temp/{TOKEN}/CNFsearch*_complete/* apps/RaptorX/;   
         rm -r data/temp/{TOKEN}/CNFsearch*_complete/;
         export OLDWD=$PWD; cd apps/RaptorX/; ./setup.pl; cd $OLDWD;
         """

rule init_modeller:
    input:
         config["license_key_in"]
    output:
         "apps/RaptorX/modeller/config.py"
    priority:
         50
    shell:
         """
         export key=$(head -n 1 {input});
         export config_file=$(mod9.23 -v 2> >(grep config.py));
         perl -i -pe "s/XXXX/$key/g" $config_file;
         export stderr_out=$(mod9.23 -v 2> >(grep 'Cannot open file -v')); 
         rm -- -v.log;
         if [ "$stderr_out" = "" ]; then 
            echo "Check your MODELLER license key.";
         else
            cp $config_file {output}
         fi
         """

rule download_databases:
    input:
         f"data/temp/{TOKEN}/{{database}}_download_link.txt"
    output:
         "apps/RaptorX/databases/archives/{database}.tar.gz"
    priority:
         40
    threads:
         1000
    shell:
         """
         wget $(head -n 1 {input}) -O {output} -q --show-progress --progress=bar:force:noscroll
         """

rule unzip_databases:
    input:
         "apps/RaptorX/databases/archives/{database}.tar.gz"
    output:
         "apps/RaptorX/databases/{database}/{database}.unzipped.txt"
    priority:
         30
    threads:
         1000
    run:
         db = wildcards.database
         if db in ["nr70", "nr90"]:
             target_dir = f"apps/RaptorX/databases/" + db
         else:
             target_dir = f"apps/RaptorX/databases/"
         shell(f"""tar -zxf {{input}} -C {target_dir};
                   touch {str(output)}""")

rule move_databases:
    input:
         "apps/RaptorX/databases/{database}/{database}.unzipped.txt"
    output:
         "apps/RaptorX/databases/{database}/{database}.moved.txt"
    priority:
         20
    run:
         db = wildcards.database
         move_files = True

         if "Remain" in db:
             target = db.replace("Remain", "BC100")
             target_dir = f"apps/RaptorX/databases/{target}/"
             shell(f"mkdir -p {target_dir}")
         elif "nr" in db:
             target_dir = "apps/RaptorX/databases/NR_new/"
             shell("mkdir -p apps/RaptorX/databases/NR_new")
         else:
             move_files = False
             print(f"Got {db}. Nothing to do.")

         if move_files:
            shell(f"""for file in `ls -1 apps/RaptorX/databases/{db}/`; do
                            mv apps/RaptorX/databases/{db}/$file {target_dir};
                      done""")

         shell("touch {output}")
