FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate plots_sequence_length_distribution" >> ~/.bashrc
ENV PATH /opt/conda/envs/plots_sequence_length_distribution/bin:\$PATH
RUN echo \$PATH
RUN conda env list
VOLUME /snakemake
WORKDIR /snakemake
ENTRYPOINT ["snakemake"]