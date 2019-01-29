FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# install CopraRNA2
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install --yes -c bioconda coprarna && \
    conda clean -ya

# update supported organisms
COPY CopraRNA/update_kegg2refseq/run/ /opt/conda/bin/coprarna_aux/

WORKDIR /mnt
CMD [ "/bin/bash" ]
