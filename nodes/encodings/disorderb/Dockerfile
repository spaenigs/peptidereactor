FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate encodings_disorderb" >> ~/.bashrc
ENV PATH /opt/conda/envs/encodings_disorderb/bin:\$PATH
RUN echo \$PATH
RUN conda env list
VOLUME /snakemake
WORKDIR /snakemake
ENTRYPOINT ["snakemake"]