FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

ENV LANG=C.UTF-8 LC_ALL=C

# install coprarna2 
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -c bioconda coprarna && \
    conda clean -ya 

# update supported organisms
COPY CopraRNA/update_kegg2refseq/run/ /opt/conda/bin/coprarna_aux/


WORKDIR /mnt
CMD [ "/bin/bash" ]
