#!/usr/bin/env bash

docker run -v /home/spaenigs/PycharmProjects/eb/nodes/utils/split_normalize/snakemake:/snakemake \
           -u `id -u $USER` \
           node_utils_split_normalize_image \
           $1