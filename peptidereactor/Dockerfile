FROM continuumio/miniconda3

RUN echo "deb http://security.debian.org/debian-security jessie/updates main" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y libgfortran3

RUN apt-get update && apt-get install -y unzip gzip
RUN apt-get update && apt-get install -y r-base r-base-dev
RUN apt-get update && apt-get install -y rsync
RUN R -e "install.packages(c('Interpol', 'yaml', 'kaos', 'seqinr', 'pROC'), repos = 'https://cran.uni-muenster.de/')"

RUN apt-get update && apt-get install -y pdb2pqr apbs

ADD iFeature/ iFeature/
RUN chmod -R 777 iFeature/data/

ADD conda/environment.yaml environment.yaml
RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate peptidereactor-env" >> ~/.bashrc

ENV PATH /opt/conda/envs/peptidereactor-env/bin:$PATH

RUN echo $PATH
RUN conda env list

RUN chmod -R 777 /opt/conda/envs/peptidereactor-env/lib/modeller-9.23/modlib/modeller/

VOLUME /snakemake
WORKDIR /snakemake

ENTRYPOINT ["snakemake", "--nolock"]