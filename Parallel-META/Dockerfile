FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y libxml2-dev libsm6 libxext6 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install parallel-meta
WORKDIR /
ENV ParallelMETA=/parallel-meta
ENV PATH=$PATH:${ParallelMETA}/bin
RUN wget http://bioinfo.single-cell.cn/Released_Software/parallel-meta/3.5/parallel-meta-3.5-bin.tar.gz -O ${ParallelMETA}.tar.gz && \
    tar zxf ${ParallelMETA}.tar.gz && rm ${ParallelMETA}.tar.gz && chmod -R 777 ${ParallelMETA}

# install r via conda 
RUN conda config --add channels r && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install --yes -c r r-base && \
    conda clean -ya

# install r-packages via conda 
RUN conda install --quiet --yes \
    'r-devtools=1.13*' \
    'r-tidyverse=1.2*' \
    'r-gplots' \
    'r-vegan' \
    'r-igraph' \
    'r-reshape' \
    'r-pheatmap' && \
    conda clean -tipsy

# install r packages in r
RUN R --vanilla -e 'install.packages("getopt", repos="http://cran.r-project.org", dependences=TRUE)'
RUN R --vanilla -e 'install.packages("https://cran.r-project.org/src/contrib/optparse_1.6.1.tar.gz", repos = NULL, type = "source", clean = TRUE)'
RUN R --vanilla -e 'install.packages(c("pROC","combinat","plyr","RColorBrewer","lattice","squash","fossil","abind", "randomForest", "psych","ade4"), repos="http://cran.r-project.org", dependences=TRUE, clean=TRUE)'

# add usearch to PATH 
ADD usearch /usr/local/bin/

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

