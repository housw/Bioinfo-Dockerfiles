FROM continuumio/miniconda3:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make cmake build-essential libtool zlib1g-dev git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install fastpath
RUN git clone https://github.com/deprekate/fastpath.git && \
    chmod -R 777 fastpath && cd fastpath && make && cd ../
RUN git clone --recursive https://github.com/deprekate/PHANOTATE.git && \
    chmod -R 777 PHANOTATE && cd PHANOTATE && make && cd ../
ENV PATH=$PATH:/fastpath:PHANOTATE/

# install tRNAscan
RUN conda install -c bioconda trnascan-se && \
    conda clean -ya
#../phanotate.py NC_000866.1.fasta --outfmt genbank --outfile test.gb

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]


