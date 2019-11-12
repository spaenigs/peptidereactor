FROM continuumio/miniconda3
RUN echo "deb http://security.debian.org/debian-security jessie/updates main" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y libgfortran3 \
    && apt-get install -y unzip
ADD iFeature/ iFeature/
RUN chmod -R 777 iFeature/data/
ADD environment.yaml environment.yaml
RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate encoding_benchmark" >> ~/.bashrc
ENV PATH /opt/conda/envs/encoding_benchmark/bin:$PATH
RUN echo $PATH
RUN conda env list
VOLUME /snakemake
WORKDIR /snakemake
ENTRYPOINT ["snakemake"]