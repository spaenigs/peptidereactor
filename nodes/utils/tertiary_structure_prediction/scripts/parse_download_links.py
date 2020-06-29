import re

with open(snakemake.input[0]) as index_html:
    for line in index_html.readlines():
        if re.match(f".*?({snakemake.wildcards.database}).*", line) is not None:
            download_link = re.match(".*href='(.*?)'", line).group(1)

with open(snakemake.output[0], mode="w") as file_out:
    file_out.write(download_link)
    file_out.flush()
