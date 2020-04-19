FROM continuumio/miniconda2:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool zlib1g-dev wget zip unzip parallel && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install perl packages
RUN conda install --yes -c bioconda samtools perl-bio-samtools perl-text-nsp && \
    conda clean -ya

# install diversitools
#RUN git clone https://github.com/josephhughes/DiversiTools.git && \
#    chmod -R 777 /DiversiTools  
RUN wget https://github.com/josephhughes/DiversiTools/archive/master.zip && \
    unzip master.zip && mv DiversiTools-master DiversiTools && \
    chmod -R 777 /DiversiTools && rm master.zip
ENV PATH=$PATH:/DiversiTools/bin

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
