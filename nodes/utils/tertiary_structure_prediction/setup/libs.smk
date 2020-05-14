TOKEN = config["token"]

rule download_libs:
    input:
         "nodes/utils/tertiary_structure_prediction/license/spineXpublic.txt"
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