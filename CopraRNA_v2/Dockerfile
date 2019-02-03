FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# install CopraRNA2
ADD JensCopraRNA2.yml /tmp/

# create coprarna2 env
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda env create --name coprarna2 --file /tmp/JensCopraRNA2.yml && \
    conda clean -ya

# update default conda env and path
ENV CONDA_DEFAULT_ENV=coprarna2
ENV COPRARNA2_PREFIX=/opt/conda/envs/$CONDA_DEFAULT_ENV
ENV PATH=$COPRARNA2_PREFIX/bin:$PATH

# update supported organisms
COPY CopraRNA/update_kegg2refseq/run/ $COPRARNA2_PREFIX/bin/coprarna_aux/

WORKDIR /mnt
CMD [ "/bin/bash" ]
