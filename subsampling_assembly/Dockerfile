FROM continuumio/miniconda3
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive

# install required tools
RUN conda install --yes python=3.5 && \
    conda install --yes -c bioconda seqtk spades amos cd-hit && \
    ln -s /opt/conda/bin/show-coords /usr/local/bin/show-coords && \
    ln -s /opt/conda/bin/blat /usr/local/bin/blat
RUN conda install --yes -c conda-forge julia
ENV JULIA_DEPOT_PATH /opt/julia_lib
RUN julia -e 'using Pkg; Pkg.add("ArgParse"); using ArgParse' && \
    julia -e 'using Pkg; Pkg.add("FastaIO"); using FastaIO'

# files
COPY filter_fasta_by_length.jl /usr/bin/
COPY opt /opt/

# env
ENV PATH /opt/conda/bin:/opt/454/bin:$PATH

# permisson
RUN chmod -R 777 /opt/*

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
