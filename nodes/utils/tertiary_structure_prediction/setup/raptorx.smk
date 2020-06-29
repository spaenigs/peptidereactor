from glob import glob

import os

TOKEN = config["token"]

PATH_TO_RAPTORX_LICENSE = \
    "nodes/utils/tertiary_structure_prediction/license/raptorx.txt"
PATH_TO_MODELLER_LICENSE = \
    "nodes/utils/tertiary_structure_prediction/license/modeller.txt"

if not os.path.exists("peptidereactor/RaptorX/setup.pl") and \
        not os.path.exists("peptidereactor/RaptorX/modeller/config.py") and \
        len(glob("peptidereactor/RaptorX/databases/*/*.moved.txt")) != 6:

    rule download_index_html:
        input:
             PATH_TO_RAPTORX_LICENSE
        output:
             temp(f"data/temp/{TOKEN}/index.html")
        priority:
             1000
        script:
             "../scripts/download_index_html.py"

    rule parse_download_links:
        input:
             f"data/temp/{TOKEN}/index.html"
        output:
             temp(f"data/temp/{TOKEN}/{{database}}_download_link.txt")
        priority:
             1000
        script:
             "../scripts/parse_download_links.py"

    rule init_raptorx:
        input:
             f"data/temp/{TOKEN}/CNFsearch_download_link.txt"
        output:
             "peptidereactor/RaptorX/setup.pl",
             temp(f"data/temp/{TOKEN}/CNFsearch.zip")
        priority:
             1000
        shell:
             f"""
             wget $(head -n 1 {{input}}) \
                -O data/temp/{TOKEN}/CNFsearch.zip \
                -q --show-progress --progress=bar:force:noscroll;   
             unzip -q data/temp/{TOKEN}/CNFsearch.zip -d data/temp/{TOKEN}/; 
             mv data/temp/{TOKEN}/CNFsearch*_complete/* peptidereactor/RaptorX/;   
             rm -r data/temp/{TOKEN}/CNFsearch*_complete/;
             export OLDWD=$PWD; cd peptidereactor/RaptorX/; ./setup.pl; cd $OLDWD;
             """

    rule init_modeller:
        input:
             PATH_TO_MODELLER_LICENSE
        output:
             "peptidereactor/RaptorX/modeller/config.py"
        priority:
             1000
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
             "peptidereactor/RaptorX/databases/archives/{database}.tar.gz"
        priority:
             1000
        threads:
             1000
        shell:
             """
             wget $(head -n 1 {input}) -O {output} -q --show-progress --progress=bar:force:noscroll
             """

    rule unzip_databases:
        input:
             ancient("peptidereactor/RaptorX/databases/archives/{database}.tar.gz")
        output:
             "peptidereactor/RaptorX/databases/{database}/{database}.unzipped.txt"
        priority:
             1000
        threads:
             1000
        script:
             "../scripts/unzip_databases.py"

    rule move_databases:
        input:
             ancient("peptidereactor/RaptorX/databases/{database}/{database}.unzipped.txt")
        output:
             "peptidereactor/RaptorX/databases/{database}/{database}.moved.txt"
        priority:
             1000
        script:
             "../scripts/move_databases.py"
