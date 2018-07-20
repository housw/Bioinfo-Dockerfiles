FROM continuumio/miniconda3:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gcc make build-essential libtool ncbi-blast+ zlib1g-dev wget git vim python python-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get clean

# install required tools
RUN conda install --yes -c bioconda hmmer prodigal pplacer mummer 

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

