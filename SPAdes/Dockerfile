FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

RUN conda install -c bioconda spades && \
    conda clean -ya

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
