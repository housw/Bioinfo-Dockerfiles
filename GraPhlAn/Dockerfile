FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# install metaphlan2
RUN conda install -c bioconda graphlan export2graphlan && \
    conda clean -ya 

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

