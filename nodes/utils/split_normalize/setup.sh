#!/usr/bin/env bash

node_group="utils"
node_name="split_normalize"
conda_env=$node_group\_$node_name
image_name=node_$conda_env\_image

touch environment.yaml
cat > environment.yaml <<- EOF
name: $conda_env
channels:
  - bioconda
  - anaconda
  - conda-forge
dependencies:
  - biopython=1.74
  - numpy=1.17.1
  - pip=19.2.3
  - python=3.7.4
  - scipy=1.3.1
  - snakemake-minimal=5.5.4
  - pip:
    - modlamp
EOF

touch Dockerfile
cat > Dockerfile <<- EOF
FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate $conda_env" >> ~/.bashrc
ENV PATH /opt/conda/envs/$conda_env/bin:\$PATH
RUN echo \$PATH
RUN conda env list
VOLUME /snakemake
WORKDIR /snakemake
ENTRYPOINT ["snakemake"]
EOF

if docker image inspect $image_name > /dev/null 2>&1; then
  echo "Image $image_name already existing. Exit..."
else
  echo "Building image $image_name..."
  docker build -t $image_name .
  mv Dockerfile docker/
  mv environment.yaml snakemake/environment.yaml
fi



