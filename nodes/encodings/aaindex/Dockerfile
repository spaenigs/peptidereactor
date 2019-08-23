FROM continuumio/miniconda3

ADD data/environment.yaml environment.yaml

RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate encode_param_free" >> ~/.bashrc

ENV PATH /opt/conda/envs/encode_param_free/bin:$PATH

RUN echo $PATH
RUN conda env list

ENTRYPOINT ["snakemake"]
