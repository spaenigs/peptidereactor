#!/usr/bin/env bash

## $1: current nodes' working dir
## $2: root directory (e.g. /home/spaenigs/PycharmProjects/eb)
## $3: path to ncbi protein database
## $4: snakemake command

docker run -v $1/snakemake:/snakemake \
           -v $2/data:/snakemake/data \
           -v $3:/snakemake/db \
           -u `id -u $USER` \
           $IMAGE_NAME \
           $4