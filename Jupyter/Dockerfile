# shengwei/datascience-notebook 
#
# VERSION       latest

FROM jupyter/datascience-notebook
MAINTAINER Shengwei Hou, housw2010@gmail.com

USER root

# system update
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -yq dist-upgrade \
 && apt-get install -yq --no-install-recommends \
    gcc \
    make \
    build-essential \
    zlib1g-dev \
    liburi-escape-xs-perl \
    cpanminus \
    ncbi-blast+ \
    hmmer \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# perl lib
RUN cpanm URI::Escape


USER $NB_UID

# update root env 
COPY packages/packages.yml .
RUN conda env update --name root --file packages.yml
RUN rm *.yml

# install julia packages
RUN julia -e 'Pkg.init()' && \
    julia -e 'Pkg.update()' && \
    julia -e 'Pkg.add("FastaIO")' && \
    julia -e 'Pkg.add("ArgParse")'
