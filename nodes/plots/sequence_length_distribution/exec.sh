#!/usr/bin/env bash

docker run -v /home/spaenigs/PycharmProjects/eb/nodes/plots/sequence_length_distribution/snakemake:/snakemake \
           -u `id -u $USER` \
           $IMAGE_NAME \
           $1