FROM continuumio/miniconda3
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive

# install required tools
RUN conda install --yes python=3.5 && \
    conda install --yes -c bioconda seqtk seqkit spades amos cd-hit && \
    conda clean -ya && \
    ln -s /opt/conda/bin/show-coords /usr/local/bin/show-coords && \
    ln -s /opt/conda/bin/blat /usr/local/bin/blat

# copy newbler executables
COPY opt /opt/

# env
ENV PATH=/opt/conda/bin:/opt/454/bin:$PATH

# permisson
RUN chmod -R 777 /opt/*

# reinstall 64-bit mummer
RUN apt-get update && apt install --yes build-essential && \
    cd /opt/conda/opt/mummer-3.23 && \
    make clean && make CPPFLAGS="-O3 -DSIXTYFOURBITS" && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# patch mummerplot
# https://github.com/bioconda/bioconda-recipes/issues/1254
RUN perl -i -pe 's/defined \(%/\(%/' /opt/conda/bin/mummerplot

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
