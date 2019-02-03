FROM continuumio/miniconda3
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive

# install required tools
#RUN conda install --yes python=3.5 && \
RUN conda install --yes -c bioconda seqtk seqkit spades amos blast-legacy && \
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

# install cd-hit with tag MAX_SEQ=10000000 to increase memory
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz && \
    tar zxfv cd-hit-v4.6.8-2017-1208-source.tar.gz && mv cd-hit-v4.6.8-2017-1208 cd_hit && \
    cd cd_hit && make MAX_SEQ=10000000 && \
    cd cd-hit-auxtools && make && \
    rm /cd-hit-v4.6.8-2017-1208-source.tar.gz && \
    chmod -R 777 /cd_hit
ENV PATH=/cd_hit:/cd_hit/cd-hit-auxtools:/cd_hit/psi-cd-hit:$PATH


# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
