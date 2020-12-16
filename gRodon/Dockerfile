FROM shengwei/ubuntu1604py37r351:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y  libgfortran3 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install devtools, biocmanager, matrixstats, biostrings, cordon
RUN conda install --yes -c conda-forge r-devtools r-biocmanager r-matrixstats && \
    conda install --yes -c bioconda bioconductor-biostrings bioconductor-cordon && \
    conda clean -ya && \
    export TAR=/bin/tar && \
    Rscript -e 'devtools::install_github("jlw-ecoevo/gRodon")'

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
