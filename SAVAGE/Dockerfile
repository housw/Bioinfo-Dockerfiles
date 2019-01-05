FROM shengwei/ubuntu1604py27
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update --fix-missing && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# add conda channels
RUN conda config --add channels r && \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

# install required tools and libraries
RUN conda install -c bioconda -y savage && \
    conda clean -ya

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

