TOKEN = config["token"]

# rule all:
#     input:
#          expand("apps/RaptorX/databases/{database}/{database}.moved.txt",
#                 database=["nr70", "nr90"]),
#          expand("apps/RaptorX/databases/{database}/{database}.moved.txt",
#                 database=["TPL_BC40", "TPL_Remain", "TemplateLists"]),
#          expand("apps/RaptorX/databases/{database}/{database}.moved.txt",
#                 database=["pdb_BC40", "pdb_Remain"]),
#          "apps/RaptorX/modeller_activated.txt"

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
         ancient(config["download_link_in"]),
         ancient("apps/RaptorX/setup.pl")
    output:
         "apps/RaptorX/databases/archives/{database}.tar.gz"
    priority:
         40
    threads:
         1
    run:
         with open(str(input[0])) as f:
             link = list(f.readlines())[0]

         db = wildcards.database

         shell(f"wget {link} -O data/temp/{TOKEN}/index_{db}.html")

         if db in ["nr70", "nr90"]:
             pattern = db
         elif db in ["TPL_BC40", "TPL_Remain", "pdb_BC40",
                     "pdb_Remain"]:  #, "TemplateLists"]:
             pattern = db + ".*"
         else:
             raise ValueError(f"Unknown database: {db}.")

         shell(f"""wget $(cat data/temp/{TOKEN}/index_{db}.html | grep {pattern}.tar.gz | \
                        perl -ne 'print "$1\n" if /(http.*?)\W>/') -O {str(output)} \
                        -q --show-progress --progress=bar:force:noscroll""")

rule unzip_databases:
    input:
         "apps/RaptorX/databases/archives/{database}.tar.gz"
    output:
         "apps/RaptorX/databases/{database}/{database}.unzipped.txt"
    priority:
         30
    run:
         db = wildcards.database
         if db in ["nr70", "nr90"]:
             target_dir = f"apps/RaptorX/databases/" + db
         else:
             target_dir = f"apps/RaptorX/databases/"
         shell(f"""tar -zxf {{input}} -C {target_dir};
                   touch {str(output)}""")

rule move_dbs:
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
         # elif db == "TemplateLists":
         #     target_dir = "apps/RaptorX/databases/"
         else:
             move_files = False
             print(f"Got {db}. Nothing to do.")

         if move_files:
            shell(f"""for file in `ls -1 apps/RaptorX/databases/{db}/`; do
                            mv apps/RaptorX/databases/{db}/$file {target_dir};
                      done""")

         shell("touch {output}")