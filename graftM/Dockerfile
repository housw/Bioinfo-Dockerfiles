FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y hmmer mafft fasttree && \
    ln -s /usr/bin/fasttree /usr/bin/FastTree && \
    ln -s /usr/bin/fasttreeMP /usr/bin/FastTreeMP && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install graftm
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -c bioconda diamond pplacer krona orfm graftm && \
    conda clean -ya

# install fxtract
WORKDIR /
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y git && \
    git clone --recursive https://github.com/ctSkennerton/fxtract.git && \
    cd fxtract && \
    make NO_PCRE=1 && \
    apt-get remove -y --auto-remove git && \
    apt-get clean && rm -rf /var/lib/apt/lists/* 
ENV PATH=$PATH:/fxtract/

WORKDIR /mnt
CMD [ "/bin/bash" ]

