from snakemake import shell

TOKEN = snakemake.config["token"]

with open(snakemake.input[0]) as file_in:
    link_to_index_html = list(file_in.readlines())[0].rstrip()
    shell(f"wget {link_to_index_html} -P data/temp/{TOKEN}/ 2> /dev/null")
