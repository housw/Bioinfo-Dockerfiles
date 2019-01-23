FROM continuumio/miniconda2:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool zlib1g-dev samtools && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install pplacer
RUN conda install --yes -c shengwei abawaca && \
    conda install -c bioconda perl-bioperl && \
    conda clean -ya

# add to 
ADD esom/prepare_esom_files.pl /opt/conda/bin/

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]


