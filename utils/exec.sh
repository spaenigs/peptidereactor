#!/usr/bin/env bash

### $1: current nodes' working dir
### $2: root directory (e.g. /home/spaenigs/PycharmProjects/eb)
### $3: path to ncbi protein database
### $4: snakemake command
#
#docker run -v $1/snakemake:/snakemake \
#           -v $2/data:/snakemake/data \
#           -v $3:/snakemake/db \
#           -u `id -u $USER` \
#           $IMAGE_NAME \
#           $4

docker run -v $PWD:/snakemake \
           -v $PWD/data:/snakemake/data \
           -v /media/spaenigs/4B1DB7375F3291A1/uniprot_db/uniref90:/snakemake/db \
           -v $PWD/apps/encoder:/snakemake/encoder \
           -v $PWD/apps/iFeature:/snakemake/iFeature \
           -v $PWD/apps/Interpol:/snakemake/Interpol \
           -v $PWD/apps/VSL2:/snakemake/VSL2 \
           -v $PWD/apps/spineXpublic:/snakemake/spineXpublic \
           -u `id -u $USER` \
           encoding_benchmark \
           $1

# clean up...
rm -rf db/ encoder/ iFeature/ Interpol/ VSL2/ spineXpublic/
