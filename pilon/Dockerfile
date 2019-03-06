FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get update && \
    apt-get install -yq --no-install-recommends zlib1g-dev libbz2-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install pilon
RUN conda install -c bioconda bwa samtools pilon && \
    conda clean -ya

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
