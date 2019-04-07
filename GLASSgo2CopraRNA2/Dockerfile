FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y ncbi-blast+ clustalo bioperl software-properties-common \
            r-base liblapack-dev libopenblas-dev gcc g++ gfortran && \
    add-apt-repository ppa:j-4/vienna-rna && apt-get -qq update && apt-get install -y vienna-rna && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install python package
RUN pip install numpy biopython ete3 csvtotable

# add scripts
RUN mkdir /GLASSgo && chmod -R 777 /GLASSgo
ADD scripts /GLASSgo/
ENV PATH $PATH:/GLASSgo/
WORKDIR /GLASSgo/

# update taxa.sqlite, install R packages
RUN ./createJSON-1_0.py -u true && rm taxdump.tar.gz && \
    Rscript install_dependences.R

WORKDIR /mnt
CMD [ "/bin/bash" ]

