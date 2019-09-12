#!/usr/bin/env bash

## $1: name of node group (e.g. 'utils')
## $2: node name
## $3: conda dependencies and their channels

node_group=$1
node_name=$2
conda_env=$node_group\_$node_name
image_name=node_$conda_env\_image

export IMAGE_NAME=$image_name

#touch environment.yaml
#cat > environment.yaml <<- EOF
#name: $conda_env
#$3
#EOF

#touch Dockerfile
#cat > Dockerfile <<- EOF
#FROM continuumio/miniconda3
#ADD environment.yaml environment.yaml
#RUN conda update -y conda \
#    && conda env create -f environment.yaml \
#    && echo "conda activate $conda_env" >> ~/.bashrc
#ENV PATH /opt/conda/envs/$conda_env/bin:\$PATH
#RUN echo \$PATH
#RUN conda env list
#VOLUME /snakemake
#WORKDIR /snakemake
#ENTRYPOINT ["snakemake"]
#EOF

if docker image inspect $image_name > /dev/null 2>&1; then
  echo "Image $image_name already existing. Skipping..."
else
  echo "Building image $image_name..."
  docker build -t $image_name . --force-rm=true
fi

#rm Dockerfile environment.yaml

