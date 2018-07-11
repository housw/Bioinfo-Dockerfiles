# dRep
#
# VERSION       2.2.2

FROM continuumio/miniconda3:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gcc make build-essential libtool ncbi-blast+ zlib1g-dev wget git vim python python-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get clean

# install required tools
RUN conda install --yes -c bioconda mummer prodigal hmmer pplacer 

# install dRep from PyPI
RUN pip install drep

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
RUN make
RUN make install prefix=/usr/local

# install mash v2.0
RUN wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar && \
    tar xfv mash-Linux64-v2.0.tar && \
    mv mash-Linux64-v2.0/mash /usr/local/bin/ && \
    rm -rf mash-Linux64-v2.0  mash-Linux64-v2.0.tar

# install checkM
ENV checkM_DIR /checkM
RUN mkdir -p $checkM_DIR
WORKDIR $checkM_DIR
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && /usr/bin/python get-pip.py 
RUN /usr/local/bin/pip install numpy scipy matplotlib pysam dendropy 
RUN git clone https://github.com/Ecogenomics/CheckM.git
WORKDIR $checkM_DIR/CheckM
RUN /usr/bin/python setup.py install
# configure checkM database
ENV checkM_db /db/CheckM
RUN mkdir -p $checkM_db
RUN echo -e "cat << EOF\n/db/CheckM\nEOF\n" | checkm data setRoot
WORKDIR $checkM_db
RUN wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
RUN tar zxfv checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

