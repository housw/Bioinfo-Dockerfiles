# dRep:2.2.3
FROM continuumio/miniconda3:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool zlib1g-dev python python-dev python-pip python-wheel python-setuptools ncbi-blast+ hmmer prodigal mummer wget git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install pplacer
RUN conda install --yes -c bioconda pplacer && \
    conda clean -ya

# install python packages
RUN /usr/bin/pip install numpy scipy matplotlib pysam dendropy checkm-genome
RUN /opt/conda/bin/pip install drep

# configure checkM database
ENV checkM_db /db/CheckM
RUN mkdir -p $checkM_db
#RUN echo -e "cat << EOF\n/db/CheckM\nEOF\n" | checkm data setRoot
RUN echo -e "\n/db/CheckM\n" | checkm data setRoot
WORKDIR $checkM_db
RUN wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar zxfv checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz


# install gANI v1
ENV gANI_DIR /gANI
RUN mkdir $gANI_DIR
WORKDIR $gANI_DIR
COPY ANIcalculator_v1.tgz .
RUN tar zxfv ANIcalculator_v1.tgz && \
    ln -s ANIcalculator_v1 current && \
    rm ANIcalculator_v1.tgz 
ENV PATH="$PATH":$gANI_DIR/current

# variables
ENV CENTRIFUGE_DIR /centrifuge
RUN mkdir $CENTRIFUGE_DIR
WORKDIR $CENTRIFUGE_DIR
# install centrifuge from github
RUN git clone https://github.com/infphilo/centrifuge
WORKDIR $CENTRIFUGE_DIR/centrifuge
# fix gcc flags
# https://github.com/infphilo/centrifuge/issues/108 
RUN sed -i 's/^RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3$/RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3 -std=c++14 -fpermissive/' Makefile
RUN make && make install prefix=/usr/local

# install mash v2.0
RUN wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar && \
    tar xfv mash-Linux64-v2.0.tar && \
    mv mash-Linux64-v2.0/mash /usr/local/bin/ && \
    rm -rf mash-Linux64-v2.0  mash-Linux64-v2.0.tar


# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

