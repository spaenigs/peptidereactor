## sudo docker rm $(sudo docker ps -a -q) && sudo docker rmi $(sudo docker images -q)

FROM continuumio/miniconda3

ADD data/environment.yaml environment.yaml
#ADD data/setup.sh setup.sh

#COPY data/setup.sh /usr/local/bin/

#RUN chmod +x /usr/local/bin/setup.sh

#RUN apt-get update \
#    && apt-get install -y build-essential \
#    && apt-get install -y wget

#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#    && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b \
#    && rm Miniconda3-latest-Linux-x86_64.sh

#ENV PATH=/miniconda/bin:${PATH}

RUN conda update -y conda \
    && conda env create -f environment.yaml \
    && echo "conda activate encoding_benchmark" >> ~/.bashrc

ENV PATH /opt/conda/envs/encoding_benchmark/bin:$PATH

RUN echo $PATH
RUN conda env list

ENTRYPOINT ["snakemake"]
