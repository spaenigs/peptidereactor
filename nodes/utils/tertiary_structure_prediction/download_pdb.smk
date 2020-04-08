import requests
import re


def get_ids():
    url = "https://www.rcsb.org/pdb/rest/getCurrent"
    response = requests.get(url)
    ids_handle = response.content.decode()
    matches = \
        [re.search("<PDB structureId=\"(\w+)\" />", line)
         for line in ids_handle.split("\n") if "structureId" in line]
    ids = [m.group(1).lower() for m in matches]
    return ids[:100]


rule all:
    input:
        "peptidereactor/db/cifs/downloaded.txt"

rule download:
    output:
        "peptidereactor/db/cifs/{id}.cif"
    shell:
        """
        set +e # escape strict mode to handle rsync error
        id={wildcards.id};
        folder_id="${{id:1:2}}";
        target_dir=$(dirname {output});
        rsync -rlpt -v -z -q --delete --port=33444 \
            rsync.rcsb.org::ftp_data/structures/divided/mmCIF/$folder_id/{wildcards.id}.cif.gz \
            $target_dir 2> /dev/null;
        if [ "$?" -eq "0" ]; then
            gunzip {output}.gz
        else
            echo "########## {wildcards.id} ####### $? "
            touch {output} 
        fi
        """

rule collect:
    input:
        expand("peptidereactor/db/cifs/{id}.cif", id=get_ids())
    output:
        "peptidereactor/db/cifs/downloaded.txt"
    shell:
        "touch {output}"

