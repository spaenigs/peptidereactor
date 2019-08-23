FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda update -y conda     && conda env create -f environment.yaml     && echo "conda activate utils_split_normalize" >> ~/.bashrc
ENV PATH /opt/conda/envs/utils_split_normalize/bin:$PATH
RUN echo $PATH
RUN conda env list
VOLUME /snakemake
WORKDIR /snakemake
ENTRYPOINT ["snakemake"]
