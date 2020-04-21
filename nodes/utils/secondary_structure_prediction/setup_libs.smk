import os

TOKEN = config["token"]
LINKS_DIR = os.path.dirname(config["VSL2_download_link_in"]) + "/"

def download_input(wildcards):
    if wildcards.target == "psipred":
        return config["psipred_download_link_in"]
    elif wildcards.target == "spineXpublic":
        return config["spineXpublic_download_link_in"],
    elif wildcards.target == "VSL2":
        return config["VSL2_download_link_in"]
    else:
        raise ValueError(f"Unknown lib: {wildcards.target}")

rule download_libs:
    input:
         download_input
    output:
         temp(f"data/temp/{TOKEN}/{{target}}.tar.gz")
    priority:
         50
    threads:
         1000
    shell:
         """
         wget $(head -n 1 {input[0]}) \
            -O {output} \
            -q --show-progress --progress=bar:force:noscroll;
         """

rule unzip_libs:
    input:
         f"data/temp/{TOKEN}/{{target}}.tar.gz"
    output:
         "peptidereactor/{target}/{target}_moved.txt"
    shell:
         """
         tar xzf {input} -C peptidereactor/;
         touch {output};
         """