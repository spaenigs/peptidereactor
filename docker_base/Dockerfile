## docker rm $(docker ps -a -q) && docker rmi $(docker images -q)

FROM continuumio/miniconda3

ADD data/environment.yaml environment.yaml

RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate eb_base" >> ~/.bashrc

ENV PATH /opt/conda/envs/eb_base/bin:$PATH

RUN echo $PATH
RUN conda env list

ENTRYPOINT ["snakemake"]
