FROM continuumio/miniconda2:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool zlib1g-dev hmmer wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install ruby
RUN conda install --yes -c bioconda parallel && \
    conda install ruby && \
    conda clean -ya

# install kofamscan 
ENV kofamscan_dir '/kofamscan'
RUN mkdir -p ${kofamscan_dir}/bin && \
    cd ${kofamscan_dir}/bin && \
    wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan.tar.gz && \
    tar xvzf kofamscan.tar.gz && rm kofamscan.tar.gz && \
    cd / && chmod -R 777 ${kofamscan_dir}
ENV PATH=${PATH}:/kofamscan/bin

# install python packages
#RUN pip install numpy scipy matplotlib pysam dendropy checkm-genome

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]


